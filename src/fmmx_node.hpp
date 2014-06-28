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
#include "expansion.hpp"
#include "valarray.hpp"
#include "delayed_action.hpp"
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
	exafmm_kernel<P> exafmm;
	cube_poles<P> cubes;
	fmmx_node_static_data() {
		position_array.resize(pow_<3,Ndim>::value);
		near_indexes.resize(pow_<3,Ndim>::value);
		neighbor_indexes.resize(pow_<3,Ndim>::value);
		std::valarray<std::size_t> dims(Ndim,Nx);
		auto position_array0 = create_position_array(dims);
		for( dir_type<Ndim> dir; !dir.is_end(); dir++) {
			std::valarray<double> corner;
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
			position_array[dir] = create_position_array(dims) + corner;
			near_indexes[dir].resize(product(dims));
			neighbor_indexes[dir].resize(product(dims));
			for( std::size_t i1 = 0; i1 < position_array[dir].size(); i1++) {
				std::list<std::size_t> list_near;
				std::list<std::size_t> list_neighbor;
				for( std::size_t i0 = 0; i0 < position_array[dir].size(); i0++) {
					if( i1 != i0 ) {
						auto dif1 = std::size_t(abs(position_array[dir][i1] - position_array0[i0]).max() + 0.5);
						auto dif0 = std::size_t(2.0 * abs(position_array[dir][i1] - 2.0 * position_array0[i0]).max() + 0.5)/2;
						if( dif0 < BOUND_WIDTH ) {
							list_neighbor.push_back(i0);
							if( dif1 < BOUND_WIDTH ) {
								list_near.push_back(i0);
							}
						}
					}
				}
				near_indexes[i1].resize(list_near.size());
				neighbor_indexes[i1].resize(list_neighbor.size());
				std::copy(std::begin(list_near), std::end(list_near), std::begin(near_indexes[i1]));
				std::copy(std::begin(list_neighbor), std::end(list_neighbor), std::begin(neighbor_indexes[i1]));
			}
		}
	}
};

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
class fmmx_node: public node<fmmx_node<Ndim, Nx, P>, Ndim>, public hpx::components::managed_component_base<fmmx_node<Ndim, Nx, P>> {
public:
	static constexpr std::size_t Bw = BOUND_WIDTH;
	static constexpr std::size_t Nchild = 1 << Ndim;
	static constexpr std::size_t Size = pow_<Nx, Ndim>::value;
	using multipoles_type = std::valarray<expansion<P>>;
	using expansions_type = std::valarray<expansion<P>>;
	using exchange_type = delayed_action<multipoles_type>;
	using descend_type = delayed_action<multipoles_type>;
	using ascend_type = delayed_action<expansions_type>;
	using type_holder = fmmx_node;
	using base_type_holder = node<fmmx_node, Ndim>;
	using wrapped_type = fmmx_node;
private:
	static fmmx_node_static_data<Ndim, Nx, P> static_data;
	std::valarray<double> rho;
	std::valarray<expansion<P>> M;
	std::valarray<expansion<P>> L;
	std::valarray<std::valarray<double>> X;
	std::valarray<double> dx;
public:
	std::vector<ascend_type> ascend(ascend_type& parent) {
		std::valarray<std::size_t> dims(Ndim, Nx);
		auto dir_center = dir_type<Ndim>().set_zero();
		exchange_type center = exchange_get(dir_center);
		exchange_set(dir_center, center);
		auto Lcoarse = parent.get();
		decltype(Lcoarse) Lthis(Lcoarse.size());
		for (std::size_t ci = 0; ci < Nchild; ci++) {
			std::fill(std::begin(Lthis), std::end(Lthis), 0.0);
			std::valarray<double> dist(Ndim);
			for (std::size_t d = 0; d < Ndim; d++) {
				if ((ci >> d) & 1) {
					dist[ci] = +0.25 * dx[d];
				} else {
					dist[ci] = -0.25 * dx[d];
				}
			}
			std::transform(std::begin(Lcoarse), std::end(Lcoarse), std::begin(Lthis), [dist](const expansion<P>& Lin) {
				expansion<P> Lout;
				static_data.exafmm.L2L((std::valarray<complex>&)Lout, (std::valarray<complex>)Lin, dist);
				return Lout;
			});
			L[get_restrict_slice(dims, ci)] += Lthis;
		}
		std::vector<ascend_type> children;
		if (!this->is_terminal()) {
			children.resize(Nchild);
			for (std::size_t ci = 0; ci < Nchild; ci++) {
				children[ci] = ascend_type(L, [dims,ci](const expansions_type& Lc) {
					return Lc[get_xtant_slice(dims, ci)];
				});
			}
		}
		return children;

	}
	descend_type descend(std::vector<descend_type>& children) {
		std::valarray<std::size_t> dims(Ndim, Nx);
		std::fill(std::begin(L), std::end(L), 0.0);
		if (children.size() > 0) {
			for (std::size_t ci = 0; ci < Nchild; ci++) {
				M[get_xtant_slice(dims, ci)] = children[ci].get();
			}
		} else {
			auto cube = static_data.cubes.get_M(dx[0]);
			std::transform(std::begin(rho), std::end(rho), std::begin(M), [&cube](double rhoin) {
				return std::valarray<complex>(cube * complex(rhoin, 0.0));
			});
		}
		return descend_type(M, [dims,this](const multipoles_type& Mfine) {
			multipoles_type Mthis(Size/Nchild);
			multipoles_type Mcoarse(Size/Nchild);
			std::fill(std::begin(Mcoarse), std::end(Mcoarse), 0.0);
			for( std::size_t ci = 0; ci < Nchild; ci++) {
				Mthis = Mfine[get_restrict_slice(dims, ci )];
				std::valarray<double> dist(Ndim);
				for( std::size_t d = 0; d < Ndim; d++ ) {
					if( (ci >> d) & 1) {
						dist[ci] = +0.25*dx[d];
					} else {
						dist[ci] = -0.25*dx[d];
					}
				}
				for( std::size_t i = 0; i <Size/Nchild; i++ ) {
					static_data.exafmm.M2M(Mcoarse[i], Mthis[i], dist);
				}
			}
			return Mcoarse;
		});

	}
	fmmx_node() :
			node<fmmx_node, Ndim>() {
		std::valarray<std::size_t> dims(Ndim, Nx);
		X.resize(Size);
		M.resize(Size);
		L.resize(Size);
		rho.resize(Size);
		dx.resize(Ndim);
		for (std::size_t i = 0; i < Size; i++) {
			X[i].resize(Ndim);
			dx[i] = 1.0 / double(Nx) / double(1 << this->get_self().get_level());
		}
	}
	void exchange_set(dir_type<Ndim> dir, exchange_type& boundary) {
		std::valarray<double> corner0(Ndim);
		for (std::size_t d = 0; d < Ndim; d++) {
			corner0[d] = this->get_self().get_position(d);
		}
		std::valarray<std::valarray<double>> Xb = (static_data.position_array[dir] * dx) + corner0;
		dir.set_zero();
		std::valarray<std::valarray<double>> Xc = (static_data.position_array[dir] * dx) + corner0;
		std::valarray<std::valarray<std::size_t>> indexes = this->is_terminal() ? static_data.near_indexes[dir] : static_data.neighbor_indexes[dir];
		std::valarray<double> x(Xb.size());
		std::valarray<double> xdif;
		std::valarray<expansion<P>> Mb;
		Mb = boundary.get();
		for (std::size_t i = 0; i < indexes.size(); i++) {
			std::valarray<expansion<P>> l(indexes[i].size());
			std::valarray<std::valarray<double>> x = Xc[indexes[i]];
			for (std::size_t j = 0; j < indexes[i].size(); j++) {
				std::valarray<double> dist = x[j] - Xb[i];
				static_data.exafmm.M2L(l[j], Mb[i], dist);
			}
			L[indexes[i]] += l;
		}

	}
	exchange_type exchange_get(dir_type<Ndim> dir) {
		return exchange_type(M, [dir](const multipoles_type& M) {
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
			return std::valarray<expansion<P>>(M[get_slice(std::valarray<std::size_t>(Ndim,Nx), mins, maxes)]);
		});
	}
}
;

template<std::size_t Ndim, std::size_t Nx, std::size_t P>
fmmx_node_static_data<Ndim, Nx, P> fmmx_node<Ndim, Nx, P>::static_data;

}
}

#endif /* FMMX_NODE_HPP_ */
