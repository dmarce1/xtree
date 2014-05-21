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
	typedef hpx::components::managed_component_base<node<Ndim>> base_type;
	typedef indexer<int_seq_const<2, Ndim>> child_index_type;
	typedef indexer<int_seq_const<3, Ndim>, int_seq_const<-1, Ndim>> dir_type;
	typedef vector<hpx::id_type, Nchild> child_array_type;
	typedef vector<hpx::id_type, Nneighbor> neighbor_array_type;
	typedef vector<vector<hpx::id_type, Nchild>, Nneighbor> niece_array_type;
	typedef vector<std::pair<int, hpx::id_type>, Nneighbor> neighbor_notify_type;

	static bool child_is_niece_of(child_index_type ci, dir_type dir) {
		for (int i = 0; i < Ndim; i++) {
			if (dir[i] == +1) {
				if (ci[i] == 0) {
					return false;
				}
			} else if (dir[i] == -1) {
				if (ci[i] == 1) {
					return false;
				}
			}
		}
		return true;
	}

private:

	bool is_leaf;
	child_array_type children;
	neighbor_array_type neighbors;
	niece_array_type nieces;
	hpx::id_type parent;

	hpx::id_type get_niece(child_index_type ci, dir_type dir) {
		child_index_type nci;
		dir_type ndi;
		for (int i = 0; i < Ndim; i++) {
			if (dir[i] == 0) {
				nci[i] = ci[i];
				ndi[i] = dir[i];
			} else if (dir[i] == +1) {
				ndi[i] = +(nci[i] ^ 0);
				nci[i] = nci[i] ^ 1;
			} else /*if( dir[i] == -1)*/{
				ndi[i] = -(nci[i] ^ 1);
				nci[i] = nci[i] ^ 1;
			}
		}
		return nieces[ndi][nci];
	}

public:

	node(hpx::id_type _parent_id = hpx::invalid_id, neighbor_array_type _neighbors = std::vector<hpx::id_type>(Nneighbor, hpx::invalid_id)) {
		neighbors = _neighbors;
		is_leaf = true;
		for (int ni = 0; ni < Nneighbor; ni++) {
			for (int ci = 0; ci < Nchild; ci++) {
				nieces[ni][ci] = hpx::invalid_id;
			}
		}
		for (int ci = 0; ci < Nchild; ci++) {
			children[ci] = hpx::invalid_id;
		}
		parent = _parent_id;
	}

	hpx::future<void> refine() {

		is_leaf = false;

		/* Get futures for new node locations */
		auto fut0 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this]() {
			vector<hpx::future<hpx::id_type>> futures(Nchild);
			for (int ci = 0; ci < Nchild; ci++) {
				futures[ci] = server::increment_load();
			}
			return futures;
		}), std::move(hpx::make_ready_future()));

		/* Allocate new nodes on localities */
		auto fut1 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](vector<hpx::future<hpx::id_type>> ids) {
			vector<hpx::future<hpx::id_type>> futures(Nchild);
			for (child_index_type ci; !ci.end(); ci++) {
				neighbor_array_type pack;
				for( dir_type dir; !dir.end(); dir++ ) {
					pack[dir] = get_niece(ci,dir);
				}
				futures[ci] = hpx::new_<node>(ids[ci].get(), base_type::get_gid(), pack);
			}
			return futures;
		}), std::move(fut0));

		/* Assign gids to children[], notify neighbors and children */
		auto fut2 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](vector<hpx::future<hpx::id_type>> vfut) {

			for (int ci = 0; ci < Nchild; ci++) {
				children[ci] = vfut[ci].get();
			}

			vector<hpx::future<void>> futures(Nneighbor);
			for( dir_type si; !si.end(); si++) {
				if( neighbors[si] != hpx::invalid_id ) {
					futures[si] = hpx::async<action_notify_refine>(neighbors[si], si, children);
				} else {
					futures[si] = hpx::make_ready_future();
				}
			}

			return futures;

		}), std::move(fut1));

		return when_all(fut2);
	}

	void destroy_children() {

		if (!is_leaf) {
			is_leaf = true;
			for (child_index_type ci; !ci.end(); ci++) {
				children[ci] = hpx::invalid_id;
			}
		}
	}

	hpx::future<void> derefine() {
		destroy_children();
		vector<hpx::future<void>> futures(Nneighbor);
		for (dir_type dir; !dir.end(); dir++) {
			if (neighbors[dir] != hpx::invalid_id) {
				futures[dir] = hpx::async<action_notify_derefine>(neighbors[dir], dir);
			} else {
				futures[dir] = hpx::make_ready_future();
			}
		}
		return when_all(futures);
	}

	virtual ~node() {
	}

	hpx::future<void> notify_refine(dir_type dir, const child_array_type& nephews) {
		hpx::future<void> ret_fut;
		dir.flip();
		nieces[dir] = nephews;
		if (!is_leaf) {
			std::vector<hpx::future<void>> futs(Nchild / 2);
			int index = 0;
			for (child_index_type ci; !ci.end(); ci++) {
				if (child_is_niece_of(ci, dir)) {
					futs[index++] = hpx::async<action_notify_of_neighbor>(children[ci], dir, get_niece(ci, dir));
				}
			}
			futs.resize(index);
			ret_fut = when_all(futs);
		} else {
			ret_fut = hpx::make_ready_future();
		}
		return ret_fut;
	}

	hpx::future<void> notify_derefine(dir_type dir) {
		hpx::future<void> ret_fut;
		dir.flip();
		nieces[dir] = std::vector<hpx::id_type>(Nneighbor, hpx::invalid_id);
		if (!is_leaf) {
			std::vector<hpx::future<void>> futs(Nchild / 2);
			int index = 0;
			for (child_index_type ci; !ci.end(); ci++) {
				if (child_is_niece_of(ci, dir)) {
					futs[index++] = hpx::async<action_notify_of_neighbor>(children[ci], dir, hpx::invalid_id);
				}
			}
			futs.resize(index);
			ret_fut = when_all(futs);
		} else {
			ret_fut = hpx::make_ready_future();
		}
		return ret_fut;
	}

	void notify_of_neighbor(dir_type dir, hpx::id_type id) {
		neighbors[dir] = id;
	}

	using action_notify_of_neighbor = typename HPX_MAKE_ACTION(node::notify_of_neighbor)::type;
	using action_notify_refine = typename HPX_MAKE_ACTION(node::notify_refine)::type;
	using action_notify_derefine = typename HPX_MAKE_ACTION(node::notify_derefine)::type;

};
}

#endif /* NODE_HPP_ */
