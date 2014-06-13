/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include "xtree.hpp"
#include <valarray>

namespace xtree {

using iter_type = std::size_t;

template<std::size_t Ndim, std::size_t Nx>
class valnode;

template<std::size_t Ndim, std::size_t Nx, std::size_t Bw>
struct valnode_slices {
	static constexpr std::size_t Nchild = 1 << Ndim;
	static constexpr std::size_t Nsub = pow_<3, Ndim>::value;
	static constexpr std::size_t Nshift = pow_<1 + 2 * Bw, Ndim>::value;
	std::gslice octant[Nchild];
	std::gslice staggered[Nchild];
	std::gslice isub[Nsub];
	std::gslice esub[Nsub];
	std::gslice eshift[Nshift];
	valnode_slices() {
		std::valarray<std::size_t> dims(Ndim, Nx);
		std::valarray<std::size_t> strides(Ndim, Nx);
		std::valarray<std::size_t> estrides(Ndim, Nx);
		std::valarray<std::size_t> dims_over_2 = dims / std::size_t(2);
		std::valarray<std::size_t> strides_over_2;
		std::size_t start;
		strides[0] = 1;
		estrides[0] = 1;
		for (iter_type di = 1; di < Ndim; di++) {
			strides[di] = strides[di - 1] * dims[di - 1];
			estrides[di] = estrides[di - 1] * (dims[di - 1] + 2 * Bw);
		}

		for (child_index_type<Ndim> ci; !ci.end(); ci++) {
			start = 0;
			for (iter_type di = 0; di < Ndim; di++) {
				start += strides[di] * ci[di] * dims_over_2[di];
			}
			octant[ci] = std::gslice(start, dims_over_2, strides);
		}

		strides_over_2 = strides / std::size_t(2);
		for (child_index_type<Ndim> ci; !ci.end(); ci++) {
			start = 0;
			for (iter_type di = 0; di < Ndim; di++) {
				start += strides[di] * ci[di];
			}
			staggered[ci] = std::gslice(start, dims_over_2, strides_over_2);
		}

		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			start = 0;
			for (iter_type di = 0; di < Ndim; di++) {
				switch (dir[di]) {
				case 0:
					dims[di] = Nx;
					break;
				case +1:
					dims[di] = Bw;
					start += strides[di] * (dims[di] - Bw);
					break;
				case -1:
					dims[di] = Bw;
					break;
				}
			}
			isub[dir] = std::gslice(start, dims, strides);
		}

		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			start = 0;
			for (iter_type di = 0; di < Ndim; di++) {
				switch (dir[di]) {
				case 0:
					dims[di] = Nx;
					start += Bw * estrides[di];
					break;
				case +1:
					dims[di] = Bw;
					start += estrides[di] * (dims[di] + Bw);
					break;
				case -1:
					dims[di] = Bw;
					break;
				}
			}
			esub[dir] = std::gslice(start, dims, estrides);
		}

		dims = std::valarray<std::size_t>(Ndim, Nx);
		for (indexer<Ndim, 1 + 2 * Bw, -3> si; !si.end(); si++) {
			start = 0;
			for (iter_type di = 0; di < Ndim; di++) {
				start += (Bw + si[di]) * estrides[di];
			}
			eshift[si] = std::gslice(start, dims, estrides);
		}

	}
};

struct valarray_vector: public std::valarray<std::valarray<double>> {
	void resize(std::size_t Ndim, std::size_t Size) {
		this->std::valarray<std::valarray<double>>::resize(Ndim);
		for (iter_type di = 0; di < Ndim; di++) {
			(*this)[di].std::valarray<double>::resize(Size);
		}
	}
};

struct multipole {
	std::valarray<double> m;
	std::valarray<std::valarray<double>> x;
	void resize(std::size_t ndim, std::size_t sz) {
		m.resize(sz);
		x.resize(ndim);
		for (iter_type d = 0; d < ndim; d++) {
			x[d].resize(sz);
		}
	}
	template<typename Arc>
	void serialize(Arc& ar, const unsigned v) {
		boost::serialization::serialize(ar, m, v);
		for (iter_type di = 0; di < x.size(); di++) {
			boost::serialization::serialize(ar, x[di], v);
		}
	}
};

struct expansion {
	std::valarray<double> l0;
	std::valarray<std::valarray<double>> l1;
	void resize(std::size_t ndim, std::size_t sz) {
		l0.resize(sz);
		l1.resize(ndim);
		for (iter_type d = 0; d < ndim; d++) {
			l1[d].resize(sz);
		}
	}
	template<typename Arc>
	void serialize(Arc& ar, const unsigned v) {
		boost::serialization::serialize(ar, l0, v);
		for (iter_type di = 0; di < l1.size(); di++) {
			boost::serialization::serialize(ar, l1[di], v);
		}
	}
};

template<std::size_t Ndim, std::size_t Nx>
class valnode: public node<valnode<Ndim, Nx>, Ndim>, public hpx::components::managed_component_base<valnode<Ndim, Nx>> {
public:
	static constexpr std::size_t Bw = 2;
	static constexpr std::size_t Nchild = 1 << Ndim;
	static constexpr std::size_t Size = pow_<Nx, Ndim>::value;
	static constexpr std::size_t Esize = pow_<Nx + 2 * Bw, Ndim>::value;
	using type_holder = valnode;
	using base_type_holder = node<valnode, Ndim>;
	using wrapped_type = valnode;
private:
	valnode_slices<Ndim, Nx, Bw> slices;
	multipole M;
	expansion L;
	valarray_vector Xc;
	std::vector<multipole> Mb;
public:
	valnode() :
			node<valnode, Ndim>() {
		M.resize(Ndim, Size);
		L.resize(Ndim, Size);

	}
	multipole descend(const std::vector<multipole>& children) {
		multipole coarse;
		if (!this->is_terminal()) {
			for (iter_type ci = 0; ci < Nchild; ci++) {
				M.m[slices.octant[ci]] = std::move(children[ci].m);
				for (iter_type di = 0; di < Ndim; di++) {
					M.x[di][slices.octant[ci]] = std::move(children[ci].x[di]);
				}
			}
		}
		if (this->get_self().get_level() != 0) {
			coarse.resize(Ndim, Size / Nchild);
			coarse.m = 0.0;
			for (iter_type di = 0; di < Ndim; di++) {
				coarse.x[di] = 0.0;
			}
			M.x *= M.m;
			for (iter_type ci = 0; ci < Ndim; ci++) {
				coarse.m += M.m[slices.staggered[ci]];
				for (iter_type di = 0; di < Ndim; di++) {
					coarse.x[di] += M.x[di][slices.staggered[ci]];
				}
			}
			coarse.x /= coarse.m;
			for (iter_type di = 0; di < Ndim; di++) {
				Xc[di] = coarse.x[di];
			}
			M.x /= M.m;
		}
		return coarse;
	}
	std::vector<expansion> ascend(const expansion& parent) {
		multipole Mext;
		multipole Mext_sh;
		multipole this_M1;
		multipole this_M2;
		valarray_vector dX;
		std::valarray<double> rinv;
		std::vector<expansion> children;
		L.l0 = 0.0;
		for (iter_type di = 0; di < Ndim; di++) {
			L.l1[di] = 0.0;
		}
		Mext.resize(Ndim, Esize);
		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			if (!dir.is_zero()) {
				Mext.m[slices.esub[dir]] = Mb[dir].m;
				for (iter_type di = 0; di < Ndim; di++) {
					Mext.x[di][slices.esub[dir]] = Mb[dir].x[di];
				}
			} else {
				Mext.m[slices.esub[dir]] = M.m;
				for (iter_type di = 0; di < Ndim; di++) {
					Mext.x[di][slices.esub[dir]] = M.x[di];
				}
			}
		}
		dX.resize(Ndim, Size / Nchild);
		rinv.resize(Size / Nchild);
		Mext_sh.resize(Ndim, Esize);
		this_M2.resize(Ndim, Size / Nchild);
		this_M1.resize(Ndim, Size / Nchild);
		for (child_index_type<Ndim> ci0; !ci0.end(); ci0++) {
			this_M1.m = M.m[slices.staggered[ci0]];
			for (iter_type di = 0; di < Ndim; di++) {
				this_M1.x[di] = M.x[di][slices.staggered[ci0]];
			}
			for (dir_type<Ndim> dir; !dir.end(); dir++) {
				indexer<Ndim, 1 + 2 * Bw, -3> si;
				for (iter_type di = 0; di < Ndim; di++) {
					si[di] = 2 * dir[di];
				}
				Mext_sh.m = Mext.m[slices.eshift[si]];
				for (iter_type di = 0; di < Ndim; di++) {
					Mext_sh.x[di] = Mext.x[di][slices.eshift[si]];
				}
				for (child_index_type<Ndim> ci1; !ci1.end(); ci1++) {
					std::size_t mindist = INT_MAX;
					for (iter_type di = 0; di < Ndim; di++) {
						mindist = std::min(mindist, std::size_t(std::abs(ci1[di] - ci0[di] + 2 * dir[di])));
					}
					if (mindist > 1) {
						this_M2.m = Mext_sh.m[slices.staggered[ci0]];
						for (iter_type di = 0; di < Ndim; di++) {
							this_M2.x[di] = Mext_sh.x[di][slices.staggered[ci0]];
						}
						for (iter_type di = 0; di < Ndim; di++) {
							dX[di] = this_M1.x[di] - this_M2.x[di];
						}
						rinv = this_M1.x[0] * this_M1.x[0];
						for (iter_type di = 1; di < Ndim; di++) {
							rinv += this_M1.x[di] * this_M1.x[di];
						}
						rinv = 1.0 / sqrt(rinv);
						L.l0[slices.staggered[ci0]] += this_M2.m * rinv;
						for (iter_type di = 0; di < Ndim; di++) {
							L.l1[di][slices.staggered[ci0]] -= this_M2.m * pow(rinv, 3.0) * dX[di];
						}
					}
				}
			}
			L.l0[slices.staggered[ci0]] += parent.l0;
			for (iter_type di = 0; di < Ndim; di++) {
				L.l1[di][slices.staggered[ci0]] += parent.l1[di];
				L.l0[slices.staggered[ci0]] += parent.l1[di] * (this_M1.x[di] - Xc[di]);
			}
		}
		dX.resize(Ndim, 0);
		rinv.resize(0);
		Mext_sh.resize(Ndim, 0);
		this_M2.resize(Ndim, 0);
		this_M1.resize(Ndim, 0);
		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			Mb[dir].resize(Ndim, 0);
		}
		children.resize(Nchild);
		for (iter_type ci = 0; ci < Nchild; ci++) {
			children[ci].l0 = L.l0[slices.octant[ci]];
			for (iter_type di = 0; di < Ndim; di++) {
				children[ci].l1[di] = L.l1[di][slices.octant[ci]];
			}
		}

		return children;
	}

	multipole exchange_get(const dir_type<Ndim>& dir) {
		multipole mb;
		mb.m = M.m[slices.isub[dir]];
		for (iter_type di = 0; di < Ndim; di++) {
			mb.x[di] = M.x[di][slices.isub[dir]];
		}
		return mb;
	}
	void exchange_set(const dir_type<Ndim>& dir, const multipole& mb) {
		Mb[dir] = mb;
	}
};

}

using valnode_type = xtree::valnode<3,8>;
XTREE_INSTANTIATE(valnode_type, 3);

int hpx_main() {
	hpx::id_type root_gid = (hpx::new_<xtree::node_type>(hpx::find_here())).get();
	xtree::tree_type* tree_ptr = hpx::async<xtree::tree_type::action_get_this>(root_gid).get();
	valnode_type* root_node = tree_ptr->get_root().get();
	using dfunc = xtree::node_type::descend_function<xtree::multipole, &valnode_type::descend>;
	using eg_func = xtree::node_type::exchange_get_function<xtree::multipole, &valnode_type::exchange_get>;
	using es_func = xtree::node_type::exchange_set_function<xtree::multipole, &valnode_type::exchange_set>;
	using afunc = xtree::node_type::ascend_function<xtree::expansion, &valnode_type::ascend>;
	std::vector<valnode_type::operation_type> ops(3);
	auto dop = xtree::node_type::make_descend_operation<dfunc>();
	auto eop = xtree::node_type::make_exchange_operation<eg_func, es_func>();
	auto aop = xtree::node_type::make_ascend_operation<afunc>();
	ops[0] = dop;
	ops[1] = eop;
	ops[2] = aop;
	root_node->execute_operations(ops);
//root_node->compute_L();
	return hpx::finalize();
}
