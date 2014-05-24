/*
 * node.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef NODE_HPP_
#define NODE_HPP_

#include <hpx/lcos/broadcast.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/util/unwrapped.hpp>

#include <cassert>

#include "indexer.hpp"
#include "location.hpp"
#include "pow_.hpp"
#include "xtree.hpp"

namespace xtree {

template<typename Derived, typename Dims>
class tree;

template<typename Derived, typename Dims>
class node: public hpx::components::managed_component_base<node<Derived, Dims>> {
public:
	static constexpr int Ndim = Dims::dim();
	static constexpr int Nchild = 1 << Ndim;
	static constexpr int Nneighbor = pow_<3, Ndim>::get();
	static constexpr int Nniece = Nchild * Nneighbor;
private:

	typedef hpx::components::managed_component_base<node<Derived, Dims>> base_type;
	typedef indexer<int_seq_const<3, Ndim>, int_seq_const<-1, Ndim>> dir_type;
	typedef vector<hpx::id_type, Nchild> child_array_type;
	typedef vector<hpx::id_type, Nneighbor> neighbor_array_type;
	typedef vector<vector<hpx::id_type, Nchild>, Nneighbor> niece_array_type;
	typedef vector<std::pair<int, hpx::id_type>, Nneighbor> neighbor_notify_type;

protected:
	typedef indexer<int_seq_const<2, Ndim>> child_index_type;

public:
	template<typename T>
	using get_type = T (Derived::*)(location<Ndim>);

	template<typename T>
	using set_type = void (Derived::*)(location<Ndim>, T);

	enum op_type {
		ASCEND, DESCEND, EXCHANGE
	};

private:

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
	location<Ndim> self;
	boost::atomic<int> ascension_counter;
	tree<Derived, Dims>* server;

private:

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

	int future_index;
	std::vector<hpx::shared_future<void>> future_locks;

public:

	node() {
		assert(false);
	}

	node(tree<Derived, Dims>* _server, const location<Ndim>& _loc, hpx::id_type _parent_id, neighbor_array_type _neighbors) {
		server = _server;
		self = _loc;
		ascension_counter = 0;
		future_index = 1;
		future_locks.resize(1);
		future_locks[0] = hpx::make_ready_future();
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

	hpx::future<void> branch() {

		is_leaf = false;

		/* Get futures for new node locations */
		auto fut0 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this]() {
			vector<hpx::future<hpx::id_type>> futures(Nchild);
			for (int ci = 0; ci < Nchild; ci++) {
				futures[ci] = load_balancer::increment_load();
			}
			return futures;
		}), std::move(hpx::make_ready_future()));

		/* Allocate new nodes on localities */
		auto fut1 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](std::vector<hpx::future<hpx::id_type>> ids) {
			vector<hpx::future<hpx::id_type>> futures(Nchild);
			for (child_index_type ci; !ci.end(); ci++) {
				neighbor_array_type pack;
				for( dir_type dir; !dir.end(); dir++ ) {
					pack[dir] = get_niece(ci,dir);
				}
				futures[ci] = server->new_node(ids[ci].get(), self.get_child(ci), base_type::get_gid(), pack);
			}
			return futures;
		}), std::move(fut0));

		/* Assign gids to children[], notify neighbors and children */
		auto fut2 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](std::vector<hpx::future<hpx::id_type>> vfut) {

			for (int ci = 0; ci < Nchild; ci++) {
				children[ci] = vfut[ci].get();
			}

			vector<hpx::future<void>> futures(Nneighbor);
			for( dir_type si; !si.end(); si++) {
				if( neighbors[si] != hpx::invalid_id ) {
					futures[si] = hpx::async<action_notify_branch>(neighbors[si], si, children);
				} else {
					futures[si] = hpx::make_ready_future();
				}
			}

			return futures;

		}), std::move(fut1));

		return when_all(fut2);
	}

	hpx::future<void> debranch() {
		is_leaf = true;
		for (child_index_type ci; !ci.end(); ci++) {
			children[ci] = hpx::invalid_id;
		}
		std::vector<hpx::future<void>> futures(Nneighbor);
		for (dir_type dir; !dir.end(); dir++) {
			if (neighbors[dir] != hpx::invalid_id) {
				futures[dir] = hpx::async<action_notify_debranch>(neighbors[dir], dir);
			} else {
				futures[dir] = hpx::make_ready_future();
			}
		}
		return when_all(futures);
	}

	virtual ~node() {
		if (!is_leaf) {
			debranch().get();
		}
		server->delete_node(static_cast<Derived*>(this));
	}

	int get_level() const {
		return self.get_level();
	}

	hpx::future<void> notify_branch(dir_type dir, const child_array_type& nephews) {
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

	hpx::future<void> notify_debranch(dir_type dir) {
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

	Derived* get_this() {
		return static_cast<Derived*>(this);
	}


	using action_notify_of_neighbor = typename HPX_MAKE_ACTION(node::notify_of_neighbor)::type;
	using action_get_this = typename HPX_MAKE_ACTION(node::get_this)::type;
	using action_notify_branch = typename HPX_MAKE_ACTION(node::notify_branch)::type;
	using action_notify_debranch = typename HPX_MAKE_ACTION(node::notify_debranch)::type;


	template<op_type Op, typename T, get_type<T> Get, set_type<T> Set>
	void setup_op_dataflow() {
		hpx::shared_future<void> future;
		std::vector<hpx::shared_future<void>> futures;
		switch (Op) {
		case ASCEND:
			if (parent != hpx::invalid_id) {
				auto fut = hpx::async<typename operations<Op, T, Get, Set>::action_get_op_data>(parent, future_index, self);
				future = fut.then(hpx::util::unwrapped([this](T arg) {
					(static_cast<Derived*>(this)->*Set)(self.get_parent(),arg);
				}));
			} else {
				future = hpx::make_ready_future();
			}
			break;
		case DESCEND:
			if (!is_leaf) {
				futures.resize(Nchild);
				for (child_index_type ci; !ci.end(); ci++) {
					auto fut = hpx::async<typename operations<Op, T, Get, Set>::action_get_op_data>(children[ci], future_index - 1, self);
					futures[ci] = fut.then(hpx::util::unwrapped([this,ci](T arg) {
						(static_cast<Derived*>(this)->*Set)(self.get_child(ci),arg);
					}));
				}
			} else {
				futures.resize(1);
				futures[0] = hpx::make_ready_future();
			}
			future = when_all(futures).share();
			break;
		case EXCHANGE:
			futures.resize(Nneighbor);
			for (dir_type dir; !dir.end(); dir++) {
				if (neighbors[dir] != hpx::invalid_id) {
					auto fut = hpx::async<typename operations<Op, T, Get, Set>::action_get_op_data>(children[dir], future_index, self);
					futures[dir] = fut.then(hpx::util::unwrapped([this,dir](T arg) {
						(static_cast<Derived*>(this)->*Set)(self.get_neighbor(dir),arg);
					}));
				} else {
					futures[dir] = hpx::make_ready_future();
				}
			}
			future = when_all(futures).share();
			break;
		}
		future_locks[future_index] = when_all(future, future_locks[future_index - 1]).share();
		future_index++;
	}

	template<typename T, get_type<T> Get>
	hpx::future<T> get_op_data(int findex, location<Ndim> requester) {
		return future_locks[findex].then(hpx::util::unwrapped([this, requester]() {
			return (static_cast<Derived*>(this)->*Get)(requester);
		}));
	}

	template<op_type Op, typename T, get_type<T> Get, set_type<T> Set>
	struct operations {
		using action_get_op_data = typename hpx::actions::make_action<decltype(&node::get_op_data<T,Get>), &node::get_op_data<T,Get>>::type;
	};

	template <typename Arc>
	void serialize( Arc& arc, const int v ) {
		assert(false);
	}


};

}

#endif /* NODE_HPP_ */
