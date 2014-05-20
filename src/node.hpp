/*
 * node.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef NODE_HPP_
#define NODE_HPP_

#include "xtree.hpp"
#include "pow_.hpp"
#include "vector.hpp"

namespace xtree {

template<int Ndim>
class node: public hpx::components::managed_component_base<node<Ndim>> {

	static constexpr int Nchild = 1 << Ndim;
	static constexpr int Nneighbor = xtree::pow_<3, Ndim>::get();
	static constexpr int Nniece = Nchild * Nneighbor;

private:

	typedef indexer<int_seq_const<2, Ndim>> child_index_type;
	typedef indexer<int_seq_const<3, Ndim>> dir_type;
	typedef vector<hpx::id_type, Nchild> child_array_type;
	typedef vector<hpx::id_type, Nneighbor> neighbor_array_type;
	typedef vector<vector<hpx::id_type, Nchild>, Nneighbor> niece_array_type;
	typedef vector<std::pair<int, hpx::id_type>> sibling_notify_type;

	bool is_leaf;
	child_array_type children;
	neighbor_array_type neighbors;
	niece_array_type nieces;
	hpx::id_type parent;

public:

	node() {
		initialize();
	}

	void initialize() {
		is_leaf = true;
		for (int ni = 0; ni < Nniece; ni++) {
			for (int ci = 0; ci < Nchild; ci++) {
				nieces[ni][ci] = hpx::invalid_id;
			}
		}
		for (int ni = 0; ni < Nniece; ni++) {
			neighbors[ni] = hpx::invalid_id;
		}
		for (int ci = 0; ci < Nchild; ci++) {
			children[ci] = hpx::invalid_id;
		}
		parent = hpx::invalid_id;

	}

	hpx::future<void> refine() {

		auto fut1 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this]() {
			vector<hpx::future<hpx::id_type>> futures(Nchild);
			for (int ci = 0; ci < Nchild; ci++) {
				futures[ci] = hpx::new_<node>(hpx::find_here());
			}
			return futures;
		}), std::move(hpx::make_ready_future()));

		auto fut2 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](vector<hpx::future<hpx::id_type>> vfut) {
			for (int ci = 0; ci < Nchild; ci++) {
				children[ci] = vfut[ci].get();
			}
		}), std::move(fut1));

		auto fut3 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this]() {
			vector<hpx::future<void>> futures(Nneighbor);
			for( dir_type si; !si.end(); si++) {
				if( neighbors[si] != hpx::invalid_id ) {
					futures[si] = hpx::async<action_notify_refine>(neighbors[si], si, children);
				} else {
					futures[si] = hpx::make_ready_future();
				}
			}
			return futures;
		}), std::move(fut2));

		auto fut4 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](vector<hpx::future<void>> vfut) {
			for( child_index_type ci; !ci.end(); ci++) {
				sibling_notify_type pack(Nchild);
			}
		}), std::move(fut3));

		return fut4;
	}

	void derefine() {

	}

	void notify_refine(const dir_type&, const child_array_type&) {
	}
	void notify_derefine(const dir_type&) {
	}
	void notify_siblings(const sibling_notify_type&) {
	}

	using action_notify_siblings = typename HPX_MAKE_ACTION(node::notify_siblings)::type;
	using action_notify_refine = typename HPX_MAKE_ACTION(node::notify_refine)::type;
	using action_notify_derefine = typename HPX_MAKE_ACTION(node::notify_derefine)::type;

};
}

#endif /* NODE_HPP_ */
