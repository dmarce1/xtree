/*
 * fmmx_node.hpp
 *
 *  Created on: Jun 13, 2014
 *      Author: dmarce1
 */

#ifndef FMMX_NODE_HPP_
#define FMMX_NODE_HPP_

#include "exafmm.hpp"
#include "fwd.hpp"
#include "cube_poles.hpp"
#include "xtree.hpp"
#include "valarray.hpp"
#include <hpx/util/function.hpp>
#include <boost/serialization/valarray.hpp>
#include <boost/serialization/complex.hpp>

namespace xtree {
namespace fmmx {

#define BOUND_WIDTH 2

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
class fmmx_node;

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
struct fmmx_node_static_data {
	std::valarray<std::valarray<std::valarray<double>>>position_array;
	std::valarray<std::valarray<std::valarray<std::size_t>>> neighbor_indexes;
	std::valarray<std::valarray<std::valarray<std::size_t>>> near_indexes;
	std::valarray<std::valarray<std::size_t>> root_indexes;
	exafmm_kernel<P> exafmm;
	cube_poles<P> cubes;
	fmmx_node_static_data() {
		position_array.resize(pow_<3,Ndim>::value);
		near_indexes.resize(pow_<3,Ndim>::value);
		neighbor_indexes.resize(pow_<3,Ndim>::value);
		root_indexes.resize(pow_<Nx,Ndim>::value);
		std::valarray<std::size_t> dims(Nx, Ndim);
		auto position_array0 = create_position_array(dims);
		for( dir_type<Ndim> dir; !dir.is_end(); dir++) {
			position_array[dir].resize(pow_<Nx,Ndim>::value);
			std::valarray<double> corner(Ndim);
			for( std::size_t di = 0; di < Ndim; di++) {
				if( dir[di] == 0 ) {
					corner[di] = 0.0;
					dims[di] = Nx;
				} else {
					if( dir[di] > 0 ) {
						corner[di] = Nx;
					} else {
						corner[di] = -BOUND_WIDTH;
					}
					dims[di] = BOUND_WIDTH;
				}
			}
			position_array[dir] = create_position_array(dims);
			std::transform(std::begin(position_array[dir]), std::end(position_array[dir]), std::begin(position_array[dir]), [&corner](const std::valarray<double>& pin) {
						return std::valarray<double>(pin + corner);
					});
			near_indexes[dir].resize(product(dims));
			neighbor_indexes[dir].resize(product(dims));
			for( std::size_t i1 = 0; i1 < position_array[dir].size(); i1++) {
				std::list<std::size_t> list_near;
				std::list<std::size_t> list_neighbor;
				for( std::size_t i0 = 0; i0 < position_array0.size(); i0++) {
					std::valarray<int> p1(Ndim), p2(Ndim);
					for( int d = 0; d < Ndim; d++) {
						p1[d] = int(position_array[dir][i1][d]-0.5+EPS+BOUND_WIDTH);
						p2[d] = int(position_array0[i0][d]-0.5+EPS+BOUND_WIDTH);
					}
					if( std::abs(p1-p2).max() > 0 ) {
						if( std::abs(p1 / int(2) - p2 / int(2)).max() < 2 ) {
							list_near.push_back(i0);
							if( std::abs(p1-p2).max() > 1 ) {
								list_neighbor.push_back(i0);
							}
						}
					}
				}
				if( list_near.size() != 0 ) {
					near_indexes[dir][i1].resize(list_near.size());
					std::copy(std::begin(list_near), std::end(list_near), std::begin(near_indexes[dir][i1]));
					if( list_neighbor.size() != 0 ) {
						neighbor_indexes[dir][i1].resize(list_neighbor.size());
						std::copy(std::begin(list_neighbor), std::end(list_neighbor), std::begin(neighbor_indexes[dir][i1]));
					}
				}
			}
		}
		for( std::size_t i1 = 0; i1 < position_array0.size(); i1++) {
			std::list<std::size_t> list_root;
			for( std::size_t i0 = 0; i0 < position_array0.size(); i0++) {
				const auto dif = position_array0[i1]-position_array0[i0];
				const auto dif1 = std::size_t((std::abs(dif)).max() + 0.001);
				if( dif1 > 1 ) {
					list_root.push_back(i0);
				}
			}
			root_indexes[i1].resize(list_root.size());
			std::copy(std::begin(list_root), std::end(list_root), std::begin(root_indexes[i1]));
		}
	}
};

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
class fmmx_node: public node<fmmx_node<Ndim, Nx, P>, Ndim>, public hpx::components::managed_component_base<
		fmmx_node<Ndim, Nx, P>, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr> {
public:
	static constexpr std::size_t Bw = BOUND_WIDTH;
	static constexpr std::size_t Nchild = 1 << Ndim;
	static constexpr std::size_t PP = P * P;
	static constexpr std::size_t Size = pow_<Nx, Ndim>::value;
	using base_type = hpx::components::managed_component_base<fmmx_node<Ndim, Nx, P>, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr>;
	using component_type = hpx::components::managed_component<fmmx_node<Ndim, Nx, P>>;
	using multipoles_type = std::valarray<std::valarray<real>>;
	using expansions_type = std::valarray<std::valarray<real>>;
	using exchange_type = future_wrapper<multipoles_type>;
	using descend_type = future_wrapper<multipoles_type>;
	using ascend_type = future_wrapper<expansions_type>;
	using type_holder = fmmx_node;
	using base_type_holder = node<fmmx_node, Ndim>;
	using wrapped_type = fmmx_node;
private:
	static fmmx_node_static_data<Ndim, Nx, P> static_data;
	std::valarray<double> rho;
	hpx::future<void> selfi_future;
	hpx::lcos::local::spinlock lock;
	multipoles_type M;
	expansions_type L;
	std::valarray<std::valarray<double>> X;
public:
	void initialize();
	std::vector<ascend_type> ascend(ascend_type& parent);
	descend_type descend(std::vector<descend_type>& children);
	void allocate();
	fmmx_node(component_type* back_ptr);
	fmmx_node();
	virtual ~fmmx_node();bool regrid_test();
	void exchange_set(dir_type<Ndim> dir, exchange_type& boundary);
	exchange_type exchange_get(dir_type<Ndim> dir);
	std::vector<typename silo_output<Ndim>::zone> get_output_zones();
	void init_grid();
};

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
void fmmx_node<Ndim, Nx, P>::initialize() {
	std::valarray<std::size_t> dims(Nx, Ndim);
	X = create_position_array(dims);
	auto corner = this->get_self().get_position();
	for (auto i = std::begin(X); i != std::end(X); ++i) {
		for (std::size_t d = 0; d != Ndim; ++d) {
			(*i)[d] /= dims[d];
			(*i)[d] *= this->get_self().get_dx();
			(*i)[d] += corner[d];
		}
	}
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
std::vector<typename fmmx_node<Ndim, Nx, P>::ascend_type> fmmx_node<Ndim, Nx, P>::ascend(
		typename fmmx_node<Ndim, Nx, P>::ascend_type& parent) {
	printf("Ascend - %i %i %i - %i -%i\n", this->get_self().get_location(0), this->get_self().get_location(1),
			this->get_self().get_location(2), this->get_self().get_level(), hpx::get_locality_id());
	const std::valarray<double> dx(1.0 / double(Nx) / double(1 << this->get_self().get_level()), Ndim);
	std::valarray<std::size_t> dims(Nx, Ndim);
	std::valarray<real> lc(PP), lf(PP);
	expansions_type Lcoarse;
	auto Lthis = expansions_type(std::valarray<real>(Size / Nchild), P * P);
	if (this->get_level() != 0) {
		Lcoarse = parent.get();
	} else {
		Lcoarse = std::valarray<std::valarray<real>>(std::valarray<real>(Size / Nchild), P * P);
	}
	for (std::size_t ci = 0; ci != Nchild; ++ci) {
		std::valarray<double> dist(Ndim);
		for (std::size_t d = 0; d < Ndim; d++) {
			if ((ci >> d) & 1) {
				dist[d] = +0.5 * dx[d];
			} else {
				dist[d] = -0.5 * dx[d];
			}
		}
		static_data.exafmm.L2L(Lthis, Lcoarse, dist, Size / Nchild);
		boost::lock_guard<decltype(lock)> scope(lock);
		for (std::size_t p = 0; p != PP; ++p) {
			L[p][get_restrict_slice(dims, ci)] += Lthis[p];
		}
	}
	selfi_future.get();
	std::vector<typename fmmx_node<Ndim, Nx, P>::ascend_type> child_data;
	if (!this->is_terminal()) {
		child_data.resize(Nchild);
		for (std::size_t ci = 0; ci < Nchild; ci++) {
			child_data[ci] = hpx::async(hpx::launch::deferred, [dims,ci](const expansions_type& Lc) {
				auto Lf = expansions_type(std::valarray<real>(Size/Nchild), P*P);
				for(std::size_t p = 0; p != PP; ++p ) {
					Lf[p] = Lc[p][get_xtant_slice(dims, ci)];
				}
				return Lf;
			}, L);
		}
	}
	return child_data;

}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
typename fmmx_node<Ndim, Nx, P>::descend_type fmmx_node<Ndim, Nx, P>::descend(
		std::vector<typename fmmx_node<Ndim, Nx, P>::descend_type>& children) {
	printf("Descend - %i %i %i - %i\n", this->get_self().get_location(0), this->get_self().get_location(1),
			this->get_self().get_location(2), this->get_self().get_level());

	std::valarray<double> dx(1.0 / double(Nx) / double(1 << this->get_self().get_level()), Ndim);
	std::valarray<std::size_t> dims(Nx, Ndim);
	std::fill(std::begin(L), std::end(L), real(0.0));
	if (children.size() > 0) {
		for (std::size_t ci = 0; ci != Nchild; ++ci) {
			auto tmp = children[ci].get();
			for (std::size_t p = 0; p != PP; ++p) {
				M[p][get_xtant_slice(dims, ci)] = tmp[p];
			}
		}
	} else {
		auto cube = static_data.cubes.get_M(dx[0]);
		for (std::size_t p = 0; p != PP; ++p) {
			for (std::size_t i = 0; i != Size; ++i) {
				M[p][i] = cube[p] * rho[i];
			}
		}
	}
	selfi_future = hpx::async([this]() {
		typename fmmx_node<Ndim, Nx, P>::exchange_type dummy;
		auto dir_center = dir_type<Ndim>().set_zero();
		this->exchange_set(dir_center, dummy);
	});

	return hpx::async(hpx::launch::deferred, [dims,dx,this](const multipoles_type& Mfine) {
		auto Mthis = multipoles_type(std::valarray<real>(Size/Nchild), P*P);
		auto Mcoarse = multipoles_type(std::valarray<real>(Size/Nchild), P*P);
		for( std::size_t ci = 0; ci != Nchild; ++ci) {
			for( std::size_t p = 0; p!= PP; ++p) {
				Mthis[p] = Mfine[p][get_restrict_slice(dims, ci )];
			}
			std::valarray<double> dist(Ndim);
			for( std::size_t d = 0; d < Ndim; d++ ) {
				if( (ci >> d) & 1) {
					dist[d] = -0.5*dx[d];
				} else {
					dist[d] = +0.5*dx[d];
				}
			}
			constexpr auto sz = Size / Nchild;
			std::valarray<std::valarray<real>> mc(std::valarray<real>(sz), PP);
			static_data.exafmm.M2M(mc, Mthis, dist, sz);
			for( std::size_t p = 0; p != PP; ++p) {
				Mcoarse[p] += mc[p];
			}
		}
		return Mcoarse;
	}, M);

}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
void fmmx_node<Ndim, Nx, P>::allocate() {
	std::valarray<std::size_t> dims(Nx, Ndim);
	X.resize(Size);
	L = std::valarray<std::valarray<real>>(std::valarray<real>(Size), P * P);
	M = std::valarray<std::valarray<real>>(std::valarray<real>(Size), P * P);
	rho.resize(Size);
	for (std::size_t i = 0; i < Size; i++) {
		X[i].resize(Ndim);
	}
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
fmmx_node<Ndim, Nx, P>::fmmx_node(component_type* back_ptr) :
		base_type(back_ptr) {
	allocate();
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
fmmx_node<Ndim, Nx, P>::fmmx_node() :
		node<fmmx_node, Ndim>() {
	allocate();
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
fmmx_node<Ndim, Nx, P>::~fmmx_node() {
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
bool fmmx_node<Ndim, Nx, P>::regrid_test() {
	return true;
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
void fmmx_node<Ndim, Nx, P>::exchange_set(dir_type<Ndim> dir,
		typename fmmx_node<Ndim, Nx, P>::exchange_type& boundary) {
	const std::valarray<double> dx(1.0 / double(Nx) / double(1 << this->get_self().get_level()), Ndim);
	std::valarray<double> corner0(Ndim);
	for (std::size_t d = 0; d < Ndim; d++) {
		corner0[d] = this->get_self().get_position(d);
	}
	std::valarray<std::valarray<double>> Xb = (static_data.position_array[dir] * dx) + corner0;
	std::valarray<std::valarray<std::size_t>> indexes =
			(this->get_self().get_level() == 0) ?
					static_data.root_indexes :
					(this->is_terminal() ? static_data.near_indexes[dir] : static_data.neighbor_indexes[dir]);
	std::valarray<double> xdif;
	multipoles_type* Mb_ptr;
	if (!dir.is_zero()) {
		Mb_ptr = new multipoles_type(boundary.get());
		for (std::size_t i = 0; i != indexes.size(); ++i) {
			const auto sz = indexes[i].size();
			std::valarray<std::valarray<real>> l(std::valarray<real>(sz), P * P);
			std::valarray<real> mb(PP);
			std::valarray<std::valarray<real>> x = X[indexes[i]];
			std::valarray<std::valarray<real>> dist(std::valarray<real>(sz), Ndim);
			for (std::size_t d = 0; d != Ndim; ++d) {
				for (std::size_t j = 0; j != sz; ++j) {
					dist[d][j] = x[j][d] - Xb[i][d];
				}
			}
#pragma vector aligned
#pragma simd
			for (std::size_t p = 0; p != PP; ++p) {
				mb[p] = (*Mb_ptr)[p][i];
			}
			static_data.exafmm.M2L(l, mb, dist, sz);
			boost::lock_guard<decltype(lock)> scope(lock);
			for (std::size_t p = 0; p != PP; ++p) {
				L[p][indexes[i]] += l[p];
			}
		}
		delete Mb_ptr;
	} else {
		const bool is_root = this->get_self().get_level() == 0;
		auto l = static_data.exafmm.M2L_interior( M, dx[0], Nx, this->is_terminal(), is_root);
		boost::lock_guard<decltype(lock)> scope(lock);
		L += l;
	}
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
typename fmmx_node<Ndim, Nx, P>::exchange_type fmmx_node<Ndim, Nx, P>::exchange_get(dir_type<Ndim> dir) {
	return hpx::async(hpx::launch::deferred, [dir](const multipoles_type& M) {
		auto m = multipoles_type(std::valarray<real>(Size/Nchild), P*P);
		std::valarray<std::size_t> mins(Ndim), maxes(Ndim);
		for( std::size_t d = 0; d < Ndim; d++) {
			switch(dir[d]) {
				case 0:
				mins[d] = 0;
				maxes[d] = Nx;
				break;
				case +1:
				mins[d] = Nx - Bw;
				maxes[d] = Nx;
				break;
				case -1:
				mins[d] = 0;
				maxes[d] = Bw;
				break;
			}
		}
		for( std::size_t p = 0; p!= PP; ++p) {
			m[p] = M[p][get_slice(std::valarray<std::size_t>(Nx,Ndim), mins, maxes)];
		}
		return m;
	}, M);
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
std::vector<typename silo_output<Ndim>::zone> fmmx_node<Ndim, Nx, P>::get_output_zones() {
	const std::valarray<double> dx(1.0 / double(Nx) / double(1 << this->get_self().get_level()), Ndim);
	std::vector<typename silo_output<Ndim>::zone> zones(Size);
	std::valarray<std::size_t> dims(Nx, Ndim);

	for (std::size_t i = 0; i < Size; i++) {
		std::vector<double> fields(7);
		fields[0] = rho[i];
		fields[1] = L[0][i];
		fields[2] = L[1][i];
		fields[3] = L[2][i];
		fields[4] = L[3][i];
		fields[5] = M[0][i];
		fields[6] = hpx::get_locality_id();
		std::copy(std::begin(dx), std::end(dx), std::begin(zones[i].span));
		std::copy(std::begin(X[i]), std::end(X[i]), std::begin(zones[i].position));

		zones[i].fields = fields;
	}
	return zones;
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
void fmmx_node<Ndim, Nx, P>::init_grid() {
	const double R0 = 0.5;
	const std::valarray<double> x0(0.5, Ndim);
	//	printf("Initializing...\n");
	for (std::size_t i = 0; i < Size; i++) {
		auto r = std::sqrt(((X[i] - x0) * (X[i] - x0)).sum());
		if (r < R0) {
			rho[i] = 1.0;
		} else {
			rho[i] = 0.0;
		}
	}
}

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
fmmx_node_static_data<Ndim, Nx, P> fmmx_node<Ndim, Nx, P>::static_data;

}

}

#endif /* FMMX_NODE_HPP_ */
