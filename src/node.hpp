/*
 * node.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef NODE_HPP_
#define NODE_HPP_

#include <hpx/include/components.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/util/unwrapped.hpp>

#include "indexer.hpp"
#include "locality_server.hpp"
#include "location.hpp"
#include "pow_.hpp"

namespace xtree {

template<typename DerivedClass, typename Dims>
class node: public hpx::components::managed_component_base<node<DerivedClass, Dims>> {
	static constexpr int Ndim = Dims::dim();
	static constexpr int Nchild = 1 << Ndim;
	static constexpr int Nneighbor = xtree::pow_<3, Ndim>::get();
	static constexpr int Nniece = Nchild * Nneighbor;

	typedef hpx::components::managed_component_base<node<DerivedClass, Dims>> base_type;
	typedef indexer<int_seq_const<3, Ndim>, int_seq_const<-1, Ndim>> dir_type;
	typedef vector<hpx::id_type, Nchild> child_array_type;
	typedef vector<hpx::id_type, Nneighbor> neighbor_array_type;
	typedef vector<vector<hpx::id_type, Nchild>, Nneighbor> niece_array_type;
	typedef vector<std::pair<int, hpx::id_type>, Nneighbor> neighbor_notify_type;
	typedef std::set<hpx::id_type> node_dir_type;
	typedef hpx::shared_future<std::shared_ptr<void>> future_stack_type;

protected:
	typedef indexer<int_seq_const<2, Ndim>> child_index_type;

private:
	static node_dir_type leaf_node_dir;
	static node_dir_type all_node_dir;
	static hpx::id_type root_node;
	static hpx::lcos::local::mutex node_dir_mutex;

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
	std::stack<future_stack_type> future_stack;
	boost::atomic<int> ascension_counter;

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

public:

	node(const location<Ndim>& _loc = location<Ndim>(), hpx::id_type _parent_id = hpx::invalid_id,
			neighbor_array_type _neighbors = std::vector<hpx::id_type>(Nneighbor, hpx::invalid_id)) {
		self = _loc;
		ascension_counter = 0;
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

		node_dir_mutex.lock();
		leaf_node_dir.insert(this->get_gid());
		all_node_dir.insert(this->get_gid());
		if (self.get_level() == 0) {
			root_node = this->get_gid();
		}
		node_dir_mutex.unlock();
		future_stack.push(future_stack_type());

	}

	hpx::future<void> branch() {

		node_dir_mutex.lock();
		leaf_node_dir.erase(leaf_node_dir.find(this->get_gid()));
		node_dir_mutex.unlock();
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
		auto fut1 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](std::vector<hpx::future<hpx::id_type>> ids) {
			vector<hpx::future<hpx::id_type>> futures(Nchild);
			for (child_index_type ci; !ci.end(); ci++) {
				neighbor_array_type pack;
				for( dir_type dir; !dir.end(); dir++ ) {
					pack[dir] = get_niece(ci,dir);
				}
				futures[ci] = hpx::new_<node>(ids[ci].get(), self.get_child(ci), base_type::get_gid(), pack);
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
		node_dir_mutex.lock();
		leaf_node_dir.insert(this->get_gid());
		node_dir_mutex.unlock();
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
		node_dir_mutex.lock();
		all_node_dir.erase(all_node_dir.find(this->get_gid()));
		leaf_node_dir.erase(leaf_node_dir.find(this->get_gid()));
		node_dir_mutex.unlock();
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

	using action_notify_of_neighbor = typename HPX_MAKE_ACTION(node::notify_of_neighbor)::type;
	using action_notify_branch = typename HPX_MAKE_ACTION(node::notify_branch)::type;
	using action_notify_debranch = typename HPX_MAKE_ACTION(node::notify_debranch)::type;

	template<typename T, std::vector<T> (DerivedClass::*Func)(T)>
	struct action_ascend {
		typedef typename hpx::actions::make_action<decltype(&node::template ascend<T,Func>), &node::template ascend<T, Func>>::type type;
		typedef typename hpx::actions::make_action<decltype(&node::template ascend_entry<T,Func>), &node::template ascend_entry<T, Func>>::type type_entry;
		static hpx::future<void> invoke() {
			node_dir_mutex.lock();
			const int leaf_cnt = leaf_node_dir.size();
			std::vector<hpx::future<void>> futures(leaf_cnt);
			int j = 0;
			for (auto i = leaf_node_dir.begin(); i != leaf_node_dir.end(); i++) {
				child_index_type ci;
				futures[j++] = hpx::async<type_entry>(*i);
			}
			node_dir_mutex.unlock();
			return when_all(futures);
		}
	};

	template<typename T, std::vector<T> (DerivedClass::*Func)(T)>
	hpx::shared_future<void> ascend_entry() {
		hpx::shared_future<T> fut1;
		fut1 = hpx::async<typename action_ascend<T, Func>::type>(parent, self.this_child_index());
		return fut1.then(hpx::util::unwrapped([this](T args) {
			(static_cast<DerivedClass*>(this)->*Func)(args);
		}));
	}

	template<typename T, std::vector<T> (DerivedClass::*Func)(T)>
	hpx::shared_future<T> ascend(child_index_type ci) {
		hpx::future<T> fut1;
		future_stack_type& fut0 = future_stack.top();
		if (!fut0.valid()) {
			ascension_counter = Nchild;
			fut1 = hpx::async<typename action_ascend<T, Func>::type>(parent, self.this_child_index());
			fut0 = fut1.then(hpx::util::unwrapped([this](T args) {
				auto tmp0 = new std::vector<T>(std::move((static_cast<DerivedClass*>(this)->*Func)(args)));
				auto tmp1 = std::shared_ptr<std::vector<T>>(tmp0);
				return std::static_pointer_cast<void>(tmp1);
			}));
		}
		return fut0.then(hpx::util::unwrapped([this,ci,&fut0](std::shared_ptr<void> data_ptr) {
			ascension_counter--;
			if( ascension_counter == 0) {
				ascension_counter = Nchild;
				fut0 = future_stack_type();
			}
			return std::move((*std::static_pointer_cast<std::vector<T>>(data_ptr))[ci]);
		}));
		return hpx::shared_future<T>();
	}

	template<typename T, T (DerivedClass::*Func)(std::vector<T>)>
	struct action_descend {
		typedef typename hpx::actions::make_action<decltype(&node::template descend<T,Func>), &node::template descend<T, Func>>::type type;
		static hpx::future<void> invoke() {
			return hpx::async<type>(root_node);
		}
	};

	template<typename T, T (DerivedClass::*Func)(std::vector<T>)>
	hpx::shared_future<T> descend() {
		hpx::shared_future<T> rfut;
		if (!is_leaf) {
			std::vector<hpx::future<T>> child_futures;
			for (child_index_type ci; !ci.end(); ci++) {
				child_futures[ci] = hpx::async<typename action_descend<T, Func>::type>(children[ci]);
			}
			rfut = when_all(child_futures).then(hpx::util::unwrapped([this](std::vector<hpx::future<T>> child_futures) {
				std::vector<T> args;
				for( int i = 0; i < child_futures.size(); i++) {
					args[i] = child_futures[i].get();
				}
				return
				(static_cast<DerivedClass*>(this)->*Func)(args);

			}));
		} else {
			rfut = hpx::make_ready_future((static_cast<DerivedClass*>(this)->*Func)(std::vector<T>(0)));
		}
		return rfut;
	}

};

template<typename T, typename U>
hpx::id_type node<T, U>::root_node = hpx::invalid_id;

template<typename T, typename U>
typename node<T, U>::node_dir_type node<T, U>::all_node_dir;

template<typename T, typename U>
typename node<T, U>::node_dir_type node<T, U>::leaf_node_dir;

template<typename T, typename U>
hpx::lcos::local::mutex node<T, U>::node_dir_mutex;

}

#endif /* NODE_HPP_ */
