/*
 * node.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef NODE_HPP_
#define NODE_HPP_

namespace xtree {

template<typename, int, typename, op_type>
class ___setup_op_dataflow;

template<int Ndim>
class node_base {
public:
	node_base() {
	}
	virtual ~node_base() {
	}
	virtual const location<Ndim>& get_self() const = 0;
};

template<typename Derived, int Ndim>
class node: public hpx::components::abstract_managed_component_base<node<Derived, Ndim>> {
public:
	static constexpr int Nchild = pow_<2, Ndim>::value;
	static constexpr int Nneighbor = pow_<3, Ndim>::value;
	static constexpr int Nniece = Nchild * Nneighbor;
	using base_type = hpx::components::managed_component_base<Derived>;
	using child_array_type = std::array<hpx::id_type, Nchild>;
	using niece_array_type = std::array<std::array<hpx::id_type, Nchild>, Nneighbor>;
	using neighbor_notify_type = std::array<std::pair<int, hpx::id_type>, Nneighbor>;
	using neighbor_array_type = std::array<hpx::id_type, Nneighbor>;
	using wrapped_type = Derived;

private:
	static bool child_is_niece_of(child_index_type<Ndim> ci, dir_type<Ndim> dir) {
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
	hpx::id_type parent;
	location<Ndim> self;
	niece_array_type nieces;
	neighbor_array_type neighbors;
	tree<wrapped_type, Ndim>* local_tree;
	hpx::shared_future<void> last_operation_future;
	int subcycle;
	hpx::lcos::local::mutex subcycle_lock;
private:
	hpx::id_type get_niece(child_index_type<Ndim> ci, dir_type<Ndim> dir) {
		child_index_type<Ndim> nci;
		dir_type<Ndim> ndi;
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
	node() {
	}
	wrapped_type* initialize(const location<Ndim>& _loc, hpx::id_type _parent_id, neighbor_array_type _neighbors, tree<wrapped_type, Ndim>* ltree) {
		subcycle = 0;
		last_operation_future = hpx::make_ready_future();
		local_tree = ltree;
		self = _loc;
		//	cur_fut = prev_fut = hpx::make_ready_future();
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
		return dynamic_cast<wrapped_type*>(this);
	}

	virtual ~node() {
		if (!is_leaf) {
			debranch().get();
		}
		if (get_level() != 0) {
			local_tree->delete_node(static_cast<wrapped_type*>(this));
		}
	}

	template<typename Arc>
	void serialize(Arc& ar, const int) {
	}

	const location<Ndim>& get_self() const {
		return self;
	}

	bool is_terminal() const {
		return is_leaf;
	}

	hpx::future<void> branch() {
		if (!is_leaf) {
			return hpx::make_ready_future();
		}
		is_leaf = false;
		auto fut0 = hpx::make_ready_future();

		/* Allocate new nodes on localities */
		auto fut1 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](void) {
			std::vector<hpx::future<hpx::id_type>> futs(Nchild);
			for (child_index_type<Ndim> ci; !ci.end(); ci++) {
				neighbor_array_type pack;
				for( dir_type<Ndim> dir; !dir.end(); dir++ ) {
					pack[dir] = get_niece(ci,dir);
				}
				futs[ci] = local_tree->new_node(self.get_child(ci), base_type::get_gid(), pack);
			}
			return futs;
		}), std::move(fut0));

		/* Assign gids to children[], notify neighbors and children */
		auto fut2 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](std::vector<hpx::future<hpx::id_type>> vfut) {

			for (int ci = 0; ci < Nchild; ci++) {
				children[ci] = vfut[ci].get();
			}

			std::vector<hpx::future<void>> futs(Nneighbor);
			for( dir_type<Ndim> si; !si.end(); si++) {
				if( neighbors[si] != hpx::invalid_id ) {
					futs[si] = hpx::async<action_notify_branch>(neighbors[si], si, children);
				} else {
					futs[si] = hpx::make_ready_future();
				}
			}

			return futs;

		}), std::move(fut1));

		return when_all(fut2);
	}

	hpx::future<void> debranch() {
		is_leaf = true;
		std::vector<hpx::future<void>> nfutures(Nneighbor);
		std::vector<hpx::future<void>> cfutures(Nchild);
		for (child_index_type<Ndim> ci; !ci.end(); ci++) {
			if (children[ci] != hpx::invalid_id) {
				cfutures[ci] = hpx::async<action_debranch>(children[ci]);
			} else {
				cfutures[ci] = hpx::make_ready_future();
			}
		}
		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			if (neighbors[dir] != hpx::invalid_id) {
				nfutures[dir] = hpx::async<action_notify_debranch>(neighbors[dir], dir);
			} else {
				nfutures[dir] = hpx::make_ready_future();
			}
		}
		auto fut1 = when_all(cfutures).then(hpx::util::unwrapped([this](std::vector<hpx::future<void>>) {
			for (child_index_type<Ndim> ci; !ci.end(); ci++) {
				children[ci] = hpx::invalid_id;
			}
		}));
		fut1.get();
		return when_all(nfutures);
	}
	int get_level() const {
		return self.get_level();
	}

	int get_subcycle() const {
		return subcycle;
	}

	hpx::future<void> notify_branch(dir_type<Ndim> dir, const child_array_type& nephews) {
		hpx::future<void> ret_fut;
		dir.flip();
		nieces[dir] = nephews;
		if (!is_leaf) {
			std::vector<hpx::future<void>> futs(Nchild / 2);
			int index = 0;
			for (child_index_type<Ndim> ci; !ci.end(); ci++) {
				if (child_is_niece_of(ci, dir)) {
					futs[index++] = hpx::async<action_notify_of_neighbor>(children[ci], dir, get_niece(ci, dir));
				}
			}
			futs.resize(index);
			ret_fut = when_all(std::move(futs));
		} else {
			ret_fut = hpx::make_ready_future();
		}
		return ret_fut;
	}

	hpx::future<void> notify_debranch(dir_type<Ndim> dir) {
		hpx::future<void> ret_fut;
		dir.flip();
		std::fill(nieces[dir].begin(), nieces[dir].end(), hpx::invalid_id);
		if (!is_leaf) {
			std::vector<hpx::future<void>> futs(Nchild / 2);
			int index = 0;
			for (child_index_type<Ndim> ci; !ci.end(); ci++) {
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

	void notify_of_neighbor(dir_type<Ndim> dir, hpx::id_type id) {
		neighbors[dir] = id;
	}

	wrapped_type* get_this() {
		return static_cast<wrapped_type*>(this);
	}

	bound_type get_boundary_type(const dir_type<Ndim>& dir) const {
		if (neighbors[dir] != hpx::invalid_id) {
			return DECOMP;
		} else {
			for (int di = 0; di < Ndim; di++) {
				if (((dir[di] == -1) && (self.get_location(di) == 0)) || ((dir[di] == +1) && (self.get_location(di) == ((1 << self.get_level()) - 1)))) {
					return PHYS;
				}
			}
			return AMR;
		}
	}

	bool has_amr_boundary() const {
		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			if (get_boundary_type(dir) == AMR) {
				return true;
			}
		}
		return false;
	}

	void wait_my_turn(int this_subcycle) {
		subcycle_lock.lock();
		while (this_subcycle != subcycle) {
			subcycle_lock.unlock();
			hpx::this_thread::suspend();
			subcycle_lock.lock();
		}
		subcycle_lock.unlock();
	}

	template<typename T>
	using ascend_type = std::vector<T> (wrapped_type::*)(const T&);
	template<typename T>
	using descend_type = T (wrapped_type::*)(const std::vector<T>&);
	template<typename T>
	using exchange_get_type = T (wrapped_type::*)(const dir_type<Ndim>&);
	template<typename T>
	using exchange_set_type = void (wrapped_type::*)(const dir_type<Ndim>&, const T&);

	template<typename T, ascend_type<T> Function>
	struct ascend_function {
		static constexpr ascend_type<T> value = Function;
		using type = T;
	};

	template<typename T, descend_type<T> Function>
	struct descend_function {
		static constexpr descend_type<T> value = Function;
		using type = T;
	};

	template<typename T, exchange_get_type<T> Function>
	struct exchange_get_function {
		static constexpr exchange_get_type<T> value = Function;
		using type = T;
	};

	template<typename T, exchange_set_type<T> Function>
	struct exchange_set_function {
		static constexpr exchange_set_type<T> value = Function;
		using type = T;
	};

	struct operation_base {
		virtual ~operation_base() {
		}
		virtual void operator()(wrapped_type&) const = 0;
	};

	using operation_type = std::shared_ptr<operation_base>;

	void execute_operations(std::vector<operation_type> operations) {
		for (int i = 0; i < operations.size(); i++) {
			operations[i]->operator()(*this);
		}
	}

	template<typename T, typename Function>
	struct ascend_operation: operation_base {
		void operator()(wrapped_type& root) const {
			root.ascend<T, Function>(hpx::make_ready_future<T>(0), root.get_subcycle());
		}
	};

	template<typename T, typename Function>
	struct descend_operation: operation_base {
		void operator()(wrapped_type& root) const {
			root.descend<T, Function>(root.get_subcycle());
		}
	};

	template<typename T, typename Get, typename Set>
	struct exchange_operation: operation_base {
		void operator()(wrapped_type& root) const {
			root.exchange_get<T, Get, Set>(root.get_subcycle());
		}
	};

	template<typename T, typename Function>
	static operation_type make_ascend_operation() {
		auto operation = std::make_shared<ascend_operation<T, Function>>();
		return std::static_pointer_cast < operation_base > (operation);
	}

	template<typename T, typename Function>
	static operation_type make_descend_operation() {
		auto operation = std::make_shared<descend_operation<T, Function>>();
		return std::static_pointer_cast < operation_base > (operation);
	}

	template<typename T, typename Get, typename Set>
	static operation_type make_exchange_operation() {
		auto operation = std::make_shared<exchange_operation<T, Get, Set>>();
		return std::static_pointer_cast < operation_base > (operation);
	}

	template<typename T, typename Function>
	void ascend(hpx::future<T> input_data_future, int this_subcycle) {
		wait_my_turn(this_subcycle);
		auto promises = std::make_shared<std::vector<hpx::promise<T>>>();
		auto f = when_all(input_data_future, last_operation_future).then(
			hpx::util::unwrapped([this,promises](HPX_STD_TUPLE<hpx::future<T>,hpx::shared_future<void>> tupfut) {
				std::vector<T> output_data;
				auto& g = HPX_STD_GET(0, tupfut);
				auto input_data = g.get();
				output_data = (static_cast<Derived*>(this)->*(Function::value))(input_data);
				for( int i = 0; i < Nchild; i++) {
					(*promises)[i].set_value(output_data[i]);
				}
			}
		));
		if (!is_leaf) {
			for (int i = 0; i < Nchild; i++) {
				hpx::apply<action_ascend<T, Function>>(children[i], (*promises)[i].get_future(), this_subcycle);
			}
		}
		subcycle_lock.lock();
		last_operation_future = f.share();
		subcycle++;
		subcycle_lock.unlock();
	}

	template<typename T, typename Function>
	hpx::shared_future<T> descend(int this_subcycle) {
		wait_my_turn(this_subcycle);
		hpx::shared_future < T > return_future;
		hpx::future < T > f;
		if (!is_leaf) {
			std::vector < hpx::future < T >> futures(Nchild);
			for (int i = 0; i < Nchild; i++) {
				futures[i] = hpx::async<action_descend<T, Function>>(children[i], this_subcycle);
			}
			auto f = when_all(when_all(futures), last_operation_future).then(
					hpx::util::unwrapped([this](HPX_STD_TUPLE<hpx::future<std::vector<hpx::future<T>>>, hpx::shared_future<void>> tuple_futs) {
						std::vector<T> input_data;
						auto& h = HPX_STD_GET(0, tuple_futs);
						std::vector<hpx::future<T>> input_data_futures = h.get();
						input_data.resize(input_data_futures.size());
						for( int i = 0; i < input_data_futures.size(); i++ ) {
							input_data[i] = input_data_futures[i].get();
						}
						return (static_cast<Derived*>(this)->*(Function::value))(input_data);
					}));
			return_future = f.share();
		} else {
			std::vector<T> empty_vector;
			f = hpx::make_ready_future((static_cast<Derived*>(this)->*(Function::value))(empty_vector));
			return_future = f.share();
		}
		subcycle_lock.lock();
		last_operation_future = return_future;
		subcycle++;
		subcycle_lock.unlock();
		return return_future;
	}

	template<typename T, exchange_set_type<T> Set>
	void exchange_set(const dir_type<Ndim>& dir, const hpx::shared_future<T>& fut, int this_subcycle) {
		wait_my_turn(this_subcycle);
		auto f = when_all(fut, last_operation_future).then(hpx::util::unwrapped([this,dir](const T& data) {
			return (this->*Set)(dir, data);
		}));
		last_operation_future = f.share();
	}

	template<typename T, exchange_get_type<T> Get>
	void exchange_get(int this_subcycle) {
		wait_my_turn(this_subcycle);
		if (!is_leaf) {
			for (int i = 0; i < Nchild; i++) {
				hpx::apply<action_exchange_get>(children[i], this_subcycle);
			}
		}
		std::vector<hpx::shared_future<void>> futures(Nneighbor);
		for (dir_type<Ndim> dir = 0; dir < Nneighbor; dir++) {
			hpx::shared_future < T > future;
			future = last_operation_future.then(hpx::util::unwrapped([this,dir]() {
				return (this->*Get)(dir);
			}));
			if (neighbors[dir] != hpx::invalid_id) {
				futures[dir] = hpx::async<action_exchange_set>(neighbors[dir], dir, future, this_subcycle);
			} else {
				futures[dir] = hpx::make_ready_future();
			}
		}
		auto f = when_all(futures, last_operation_future);
		subcycle_lock.lock();
		last_operation_future = f.share();
		subcycle++;
		subcycle_lock.unlock();
	}

	template<typename T, typename Function>
	using action_ascend = hpx::actions::make_action< void(node::*)(hpx::future<T>, int), &node::ascend<T,Function>>; //
	template<typename T, typename Function>
	using action_descend = hpx::actions::make_action<hpx::shared_future<T>(node::*)(int), &node::descend<T,Function>>; //
	template<typename T, typename Set>
	using action_exchange_set = hpx::actions::make_action<void(node::*)(const dir_type<Ndim>&, const hpx::shared_future<T>&,int), &node::exchange_set>; //
	template<typename T, typename Get>
	using action_exchange_get = hpx::actions::make_action<void(node::*)(int), &node::exchange_get<T,Get>>; //

	HPX_DEFINE_COMPONENT_ACTION_TPL(node, initialize, action_initialize); //
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, debranch, action_debranch); //
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, get_this, action_get_this); //
	HPX_DEFINE_COMPONENT_ACTION_TPL(node,notify_branch, action_notify_branch); //
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, notify_debranch, action_notify_debranch); //
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, notify_of_neighbor, action_notify_of_neighbor);
	//
};

}

#endif /* NODE_HPP_ */
