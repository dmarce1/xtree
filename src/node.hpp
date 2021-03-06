/*
 / * node.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef NODE_HPP_
#define NODE_HPP_

#define MAX_OPS 10

namespace xtree {

template<typename, int, typename, op_type>
class ___setup_op_dataflow;

template<typename Derived, int Ndim>
class node: public hpx::components::abstract_managed_component_base<node<Derived, Ndim>> {
public:
	static constexpr int Nchild = pow_<2, Ndim>::value;
	static constexpr int Nneighbor = pow_<3, Ndim>::value;
	using base_type = hpx::components::managed_component_base<Derived>;
	using child_array_type = std::vector<hpx::id_type>;
	using niece_array_type = std::vector<std::vector<hpx::id_type>>;
	using neighbor_array_type = std::vector<hpx::id_type>;
	using wrapped_type = Derived;
	using local_type = void (wrapped_type::*)();
	using regrid_type = bool (wrapped_type::*)();
	template<typename T>
	using ascend_type = std::vector<T> (wrapped_type::*)(T&);
	template<typename T>
	using descend_type = T (wrapped_type::*)(std::vector<T>&);
	template<typename T>
	using exchange_get_type = T (wrapped_type::*)(dir_type<Ndim>);
	template<typename T>
	using exchange_set_type = void (wrapped_type::*)(dir_type<Ndim>, T&);

	template<local_type Function>
	struct local_function {
		static constexpr local_type value = Function;
		using type = void;
	};

	template<regrid_type Function>
	struct regrid_function {
		static constexpr regrid_type value = Function;
		using type = bool;
	};

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
		virtual ~operation_base() = default;
		virtual void operator()(node&) const = 0;
		virtual bool is_regrid() const;
	};
	template<typename Function>
	struct local_operation: operation_base {
		void operator()(node& root) const;
	};

	template<typename Function>
	struct regrid_operation: operation_base {
		virtual bool is_regrid() const;
		void operator()(node& root) const;
	};

	template<typename Function>
	struct ascend_operation: operation_base {
		void operator()(node& root) const;
	};

	template<typename Function>
	struct descend_operation: operation_base {
		void operator()(node& root) const;
	};

	template<typename Get, typename Set>
	struct exchange_operation: operation_base {
		void operator()(node& root) const;
	};

	using operation_type = std::shared_ptr<operation_base>;

private:
	hpx::lcos::local::counting_semaphore exchange_semaphore;
	bool is_leaf;
	bool is_branching;
	child_array_type children;
	hpx::id_type parent;
	location<Ndim> self;
	niece_array_type nieces;
	neighbor_array_type neighbors;
	tree<wrapped_type, Ndim>* local_tree;
	hpx::shared_future<void> last_operation_future;
	std::vector<hpx::promise<void>> exchange_promises;
	int subcycle;
	double myload;
	mutable hpx::lcos::local::spinlock branch_lock;
	mutable hpx::lcos::local::spinlock set_lock;
	mutable std::vector<hpx::promise<void>> subcycle_promises;
	mutable std::vector<hpx::future<void>> subcycle_futures;
private:
	hpx::id_type get_sibling_of_child(child_index_type<Ndim> ci, dir_type<Ndim> dir);
public:
	node();
	virtual ~node();
	wrapped_type* initialize(const location<Ndim>& _loc, hpx::id_type parent_id, tree<wrapped_type, Ndim>* ltree);
	template<typename Arc>
	void serialize(Arc& ar, const int);
	bool has_amr_boundary() const;
	bool is_terminal() const;
	template<typename Function>
	bool regrid(int this_subcycle);
	bound_type get_boundary_type(const dir_type<Ndim>& dir) const;
	const location<Ndim>& get_self() const;
	hpx::future<void> debranch();
	hpx::shared_future<void> operations_end(int this_subcycle);
	int get_level() const;
	int get_subcycle() const;
	wrapped_type * get_this();
	void branch();
	void execute_operations(std::vector<operation_type> operations);
	template<typename Function>
	void local(int this_subcycle);
	void wait_my_turn(int this_subcycle);
	template<typename Function>
	static operation_type make_local_operation();
	template<typename Function>
	static operation_type make_regrid_operation();
	template<typename Function>
	static operation_type make_ascend_operation();
	template<typename Function>
	static operation_type make_descend_operation();
	template<typename Get, typename Set>
	static operation_type make_exchange_operation();
	template<typename Function>
	void ascend(hpx::future<typename Function::type> input_data_future, int this_subcycle);
	template<typename Function>
	hpx::shared_future<typename Function::type> descend(int this_subcycle);
	template<typename Set>
	void exchange_set(const dir_type<Ndim>& dir, typename Set::type data);
	template<typename Get, typename Set>
	void exchange_get(int this_subcycle);
	template<typename Function>
	using action_local = hpx::actions::make_action<void(node::*)(int), &node::local<Function>>;
	template<typename Function>
	using action_regrid = hpx::actions::make_action< bool(node::*)(int), &node::regrid<Function>>;
	template<typename Function>
	using action_ascend = hpx::actions::make_action< void(node::*)(hpx::future<typename Function::type>, int), &node::ascend<Function>>;
	template<typename Function>
	using action_descend = hpx::actions::make_action<hpx::shared_future<typename Function::type>(node::*)(int), &node::descend<Function>>;
	template<typename Set>
	using action_exchange_set = hpx::actions::make_action<void(node::*)(const dir_type<Ndim>&, typename Set::type), &node::exchange_set<Set>>;
	template<typename Get, typename Set>
	using action_exchange_get = hpx::actions::make_action<void(node::*)(int), &node::exchange_get<Get,Set>>;
	void neighbors_exchange_set(dir_type<Ndim> dir, const std::vector<hpx::id_type>& nephews);
	void neighbors_ascend(std::vector<hpx::id_type>);
	void neighbors_exchange_get(int this_subcycle);
	using action_neighbors_ascend = hpx::actions::make_action<void(node::*)(std::vector<hpx::id_type>), &node::neighbors_ascend>;
	using action_neighbors_exchange_get = hpx::actions::make_action<void(node::*)(int), &node::neighbors_exchange_get>;
	using action_neighbors_exchange_set = hpx::actions::make_action<void(node::*)(dir_type<Ndim>,
			const std::vector<hpx::id_type>&), &node::neighbors_exchange_set>;
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, operations_end, action_operations_end);		//
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, initialize, action_initialize);	//
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, debranch, action_debranch);	//
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, get_this, action_get_this);
};

template<class Derived, int Ndim>
inline hpx::id_type node<Derived, Ndim>::get_sibling_of_child(child_index_type<Ndim> ci, dir_type<Ndim> dir) {
	child_index_type < Ndim > nci;
	dir_type < Ndim > ndi;
	for (int i = 0; i < Ndim; i++) {
		if (dir[i] == 0) {
			nci[i] = ci[i];
			ndi[i] = dir[i];
		} else if (dir[i] == +1) {
			ndi[i] = +(ci[i] ^ 0);
			nci[i] = ci[i] ^ 1;
		} else /*if( dir[i] == -1)*/{
			ndi[i] = -(ci[i] ^ 1);
			nci[i] = ci[i] ^ 1;
		}
	}
	return nieces[ndi][nci];
}

template<class Derived, int Ndim>
inline node<Derived, Ndim>::node() :
		exchange_semaphore(0) {
}

template<class Derived, int Ndim>
inline typename node<Derived, Ndim>::wrapped_type* node<Derived, Ndim>::initialize(const location<Ndim>& _loc,
		hpx::id_type parent_id, tree<typename node<Derived, Ndim>::wrapped_type, Ndim>* ltree) {
	subcycle_futures = std::move(std::vector<hpx::future<void>>(MAX_OPS));
	subcycle_promises.resize(MAX_OPS);
	if (parent_id != hpx::invalid_id) {
		for (int i = 0; i < MAX_OPS; i++) {
			subcycle_futures[i] = hpx::make_ready_future();
		}
	} else {
		subcycle_futures[0] = hpx::make_ready_future();
		for (int i = 1; i < MAX_OPS; i++) {
			subcycle_futures[i] = subcycle_promises[i - 1].get_future();
		}
	}
	subcycle = 0;
	last_operation_future = hpx::make_ready_future();
	local_tree = ltree;
	self = _loc;
	//	cur_fut = prev_fut = hpx::make_ready_future();
	is_leaf = true;
	is_branching = false;
	neighbors.resize(Nneighbor);
	nieces.resize(Nneighbor);
	for (int ni = 0; ni < Nneighbor; ni++) {
		nieces[ni].resize(Nchild);
		for (int ci = 0; ci < Nchild; ci++) {
			nieces[ni][ci] = hpx::invalid_id;
		}
	}
	children.resize(Nchild);
	for (int ci = 0; ci < Nchild; ci++) {
		children[ci] = hpx::invalid_id;
	}
	parent = parent_id;
	static_cast<Derived*>(this)->initialize();
	return dynamic_cast<typename node<Derived, Ndim>::wrapped_type*>(this);
}

template<class Derived, int Ndim>
inline node<Derived, Ndim>::~node() {
	if (!is_leaf) {
		debranch().get();
	}
	if (get_level() != 0) {
		local_tree->delete_node(static_cast<typename node<Derived, Ndim>::wrapped_type*>(this));
	}
	parent = hpx::invalid_id;
}

template<class Derived, int Ndim>
template<typename Arc>
inline void node<Derived, Ndim>::serialize(Arc& ar, const int) {
}

template<class Derived, int Ndim>
inline const location<Ndim>&node<Derived, Ndim>::get_self() const {
	return self;
}

template<class Derived, int Ndim>
inline bool node<Derived, Ndim>::is_terminal() const {
	return is_leaf;
}

template<class Derived, int Ndim>
inline void node<Derived, Ndim>::branch() {
	if (!is_leaf) {
		return;
	}
	is_branching = true;
	std::vector < hpx::future < hpx::id_type >> futs(Nchild);
	for (child_index_type < Ndim > ci; !ci.is_end(); ci++) {
		auto tgid = static_cast<Derived*>(this)->get_gid();
		futs[ci] = local_tree->new_node(self.get_child(ci), tgid, 0);
	}
	wait_all (futs);
	is_leaf = false;
	dir_type < Ndim > zero;
	zero.set_zero();
	for (int ci = 0; ci < Nchild; ci++) {
		children[ci] = futs[ci].get();
		nieces[zero][ci] = children[ci];
	}
//		printf("Made children\n");
}

template<class Derived, int Ndim>
inline hpx::future<void> node<Derived, Ndim>::debranch() {
	boost::lock_guard<decltype(branch_lock)> scope_lock(branch_lock);
	is_leaf = true;
	std::vector<hpx::future<void>> nfutures(Nneighbor);
	std::vector<hpx::future<void>> cfutures(Nchild);
	for (child_index_type < Ndim > ci; !ci.is_end(); ci++) {
		if (children[ci] != hpx::invalid_id) {
			cfutures[ci] = hpx::async < action_debranch > (children[ci]);
		} else {
			cfutures[ci] = hpx::make_ready_future();
		}
	}
	auto fut1 = when_all(std::move(cfutures)).then(hpx::util::unwrapped([this](std::vector<hpx::future<void>>) {
		boost::lock_guard<decltype(branch_lock)> scope_lock(branch_lock);
		for (std::size_t ci = 0; ci != Nchild; ++ci) {
			children[ci] = hpx::invalid_id;
		}
	}));
	return fut1;
}

template<class Derived, int Ndim>
inline int node<Derived, Ndim>::get_level() const {
	return self.get_level();
}

template<class Derived, int Ndim>
inline int node<Derived, Ndim>::get_subcycle() const {
	return subcycle;
}

template<class Derived, int Ndim>
inline typename node<Derived, Ndim>::wrapped_type * node<Derived, Ndim>::get_this() {
	return static_cast<typename node<Derived, Ndim>::wrapped_type*>(this);
}

template<class Derived, int Ndim>
inline bound_type node<Derived, Ndim>::get_boundary_type(const dir_type<Ndim>& dir) const {
	if (neighbors[dir] != hpx::invalid_id) {
		return DECOMP;
	} else {
		for (int di = 0; di < Ndim; di++) {
			if (((dir[di] == -1) && (self.get_location(di) == 0))
					|| ((dir[di] == +1) && (self.get_location(di) == ((1 << self.get_level()) - 1)))) {
				return PHYS;
			}
		}
		return AMR;
	}
}

template<class Derived, int Ndim>
inline bool node<Derived, Ndim>::has_amr_boundary() const {
	for (dir_type < Ndim > dir; !dir.is_end(); dir++) {
		if (get_boundary_type(dir) == AMR) {
			return true;
		}
	}
	return false;
}

template<class Derived, int Ndim>
inline void node<Derived, Ndim>::wait_my_turn(int this_subcycle) {
	assert(subcycle_futures[this_subcycle].valid());
	subcycle_futures[this_subcycle].get();
}

template<typename Derived, int Ndim>
inline bool node<Derived, Ndim>::operation_base::is_regrid() const {
	return false;
}

template<typename Derived, int Ndim>
inline hpx::shared_future<void> node<Derived, Ndim>::operations_end(int this_subcycle) {
	wait_my_turn(this_subcycle);
	hpx::future<void> future;
	if (!is_leaf) {
		std::vector<hpx::future<void>> futs(Nchild);
		for (std::size_t i = 0; i != Nchild; ++i) {
			futs[i] = hpx::async < action_operations_end > (children[i], this_subcycle);
		}
		future = when_all(std::move(futs));
	} else {
		future = hpx::make_ready_future();
	}
	auto f = when_all(future, last_operation_future).then(
			hpx::util::unwrapped([this](hpx::util::tuple<hpx::future<void>, hpx::shared_future<void>>) {
				subcycle = 0;
				subcycle_futures = std::move(std::vector<hpx::future<void>>(MAX_OPS));
				subcycle_promises.resize(MAX_OPS);
				subcycle_futures[0] = hpx::make_ready_future();
				for (int i = 1; i < MAX_OPS; i++) {
					subcycle_promises[i-1] = hpx::promise<void>();
					subcycle_futures[i] = subcycle_promises[i-1].get_future();
				}
			}));
	last_operation_future = f.share();
	//	printf("%i\n", this->get_self().get_level());

	return last_operation_future;
}

template<typename Derived, int Ndim>
inline void node<Derived, Ndim>::execute_operations(std::vector<operation_type> operations) {
	hpx::shared_future<void> fut;
	int i;
	fut = hpx::make_ready_future().share();
	for (i = 0; i < operations.size(); i++) {
		auto f = fut.then(hpx::util::unwrapped([i,this,&operations]() {
			operations[i]->operator()(*this);
		}));
		fut = f.share();
	}
	fut.get();
	if (!operations[i - 1]->is_regrid()) {
		fut = operations_end(subcycle);
	} else {
		fut = last_operation_future;
	}
	fut.get();
}

template<typename Derived, int Ndim>
template<typename Function>
inline void node<Derived, Ndim>::local_operation<Function>::operator()(node& root) const {
	root.local<Function>(root.get_subcycle());
}

template<typename Derived, int Ndim>
template<typename Function>
inline bool node<Derived, Ndim>::regrid_operation<Function>::is_regrid() const {
	return true;
}

template<typename Derived, int Ndim>
template<typename Function>
inline void node<Derived, Ndim>::regrid_operation<Function>::operator()(node& root) const {
	root.regrid<Function>(root.get_subcycle());
}

template<typename Derived, int Ndim>
template<typename Function>
inline void node<Derived, Ndim>::ascend_operation<Function>::operator()(node& root) const {
	root.ascend<Function>(hpx::make_ready_future<typename Function::type>(typename Function::type()),
			root.get_subcycle());
}

template<typename Derived, int Ndim>
template<typename Function>
inline void node<Derived, Ndim>::descend_operation<Function>::operator()(node& root) const {
	root.descend<Function>(root.get_subcycle());
}

template<typename Derived, int Ndim>
template<typename Get, typename Set>
inline void node<Derived, Ndim>::exchange_operation<Get, Set>::operator()(node& root) const {
	root.exchange_get<Get, Set>(root.get_subcycle());
}

template<typename Derived, int Ndim>
template<typename Function>
inline typename node<Derived, Ndim>::operation_type node<Derived, Ndim>::make_local_operation() {
	auto operation = std::make_shared<local_operation<Function>>();
	return std::static_pointer_cast < operation_base > (operation);
}

template<typename Derived, int Ndim>
template<typename Function>
inline typename node<Derived, Ndim>::operation_type node<Derived, Ndim>::make_regrid_operation() {
	auto operation = std::make_shared<regrid_operation<Function>>();
	return std::static_pointer_cast < operation_base > (operation);
}

template<typename Derived, int Ndim>
template<typename Function>
inline typename node<Derived, Ndim>::operation_type node<Derived, Ndim>::make_ascend_operation() {
	auto operation = std::make_shared<ascend_operation<Function>>();
	return std::static_pointer_cast < operation_base > (operation);
}

template<typename Derived, int Ndim>
template<typename Function>
inline typename node<Derived, Ndim>::operation_type node<Derived, Ndim>::make_descend_operation() {
	auto operation = std::make_shared<descend_operation<Function>>();
	return std::static_pointer_cast < operation_base > (operation);
}

template<typename Derived, int Ndim>
template<typename Get, typename Set>
inline typename node<Derived, Ndim>::operation_type node<Derived, Ndim>::make_exchange_operation() {
	auto operation = std::make_shared<exchange_operation<Get, Set>>();
	return std::static_pointer_cast < operation_base > (operation);
}

template<typename Derived, int Ndim>
template<typename Function>
inline void node<Derived, Ndim>::local(int this_subcycle) {
	wait_my_turn(this_subcycle);
	std::vector<hpx::future<void>> futures;
	hpx::shared_future<void> rc;
	if (!is_leaf) {
		futures = std::move(std::vector<hpx::future<void>>(Nchild));
		for (int i = 0; i < Nchild; ++i) {
			futures[i] = hpx::async < action_local < Function >> (children[i], this_subcycle);
		}
	}
	rc = hpx::async([this]() {
		(static_cast<Derived*>(this)->*(Function::value))();
	}).share();

	if (futures.size()) {
		rc = when_all(rc, when_all(std::move(futures))).share();
	}

	last_operation_future = rc;
	subcycle_promises[subcycle].set_value();
	subcycle++;
}

template<typename Derived, int Ndim>
template<typename Function>
inline bool node<Derived, Ndim>::regrid(int this_subcycle) {
	std::vector<hpx::future<bool>> futures(Nchild);
	bool rc;
	if (!is_leaf) {
		for (int i = 0; i < Nchild; ++i) {
			auto c = children[i];
			futures[i] = hpx::async < action_regrid < Function >> (c, this_subcycle);
		}
		wait_all(futures);
		for (auto i = 0; i != Nchild; i++) {
			if (!futures[i].get()) {
				debranch().get();
			}
		}
		rc = true;
	} else {
		rc = (static_cast<Derived*>(this)->*(Function::value))();
		if (rc) {
			branch();
		}
	}

	if (this->get_self().get_level() == 0) {
		neighbors_exchange_get(subcycle);
		printf("A\n");
		operations_end(subcycle).get();
		printf("B\n");
		neighbors_ascend(std::vector < hpx::id_type > (Nneighbor));
	}

	return rc;
}

template<typename Derived, int Ndim>
template<typename Function>
inline void node<Derived, Ndim>::ascend(hpx::future<typename Function::type> input_data_future, int this_subcycle) {
	using T = typename Function::type;
	wait_my_turn(this_subcycle);
	auto promises = std::make_shared<std::vector<hpx::promise<T>>>();
	promises->resize(Nchild);
	auto f = when_all(input_data_future, last_operation_future).then(
			hpx::util::unwrapped([=](HPX_STD_TUPLE<hpx::future<T>,hpx::shared_future<void>> tupfut) {
				std::vector<T> output_data;
				auto& g = HPX_STD_GET(0, tupfut);
				auto input_data = g.get();
				if( !is_leaf ) {
					output_data = (static_cast<Derived*>(this)->*(Function::value))(input_data);
					for( int i = 0; i < Nchild; i++) {
						(*promises)[i].set_value(output_data[i]);
					}
				} else {
					(static_cast<Derived*>(this)->*(Function::value))(input_data);
				}
			}));
	if (!is_leaf) {
		for (int i = 0; i < Nchild; i++) {
			hpx::apply < action_ascend < Function >> (children[i], (*promises)[i].get_future(), this_subcycle);
		}
	}
	last_operation_future = f.share();
	subcycle_promises[subcycle].set_value();
	subcycle++;
}

template<typename Derived, int Ndim>
template<typename Function>
inline hpx::shared_future<typename Function::type> node<Derived, Ndim>::descend(int this_subcycle) {
	//test = 0;
	using T = typename Function::type;
	wait_my_turn(this_subcycle);
	hpx::shared_future<T> return_future;
	hpx::future<T> f;
	if (!is_leaf) {
		std::vector < hpx::future < T >> futures(Nchild);
		for (int i = 0; i < Nchild; i++) {
			futures[i] = hpx::async < action_descend < Function >> (children[i], this_subcycle);
		}
		auto f =
				when_all(when_all(std::move(futures)), last_operation_future).then(
						hpx::util::unwrapped(
								[this](HPX_STD_TUPLE<hpx::future<std::vector<hpx::future<T>>>, hpx::shared_future<void>> tuple_futs) {
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

		auto f = last_operation_future.then(hpx::util::unwrapped([this]() {
			std::vector<T> empty_vector;
			return (static_cast<Derived*>(this)->*(Function::value))(empty_vector);
		}));

		return_future = f.share();
	}
	last_operation_future = return_future;
	subcycle_promises[subcycle].set_value();
	subcycle++;
	return return_future;
}

template<typename Derived, int Ndim>
template<typename Set>
inline void node<Derived, Ndim>::exchange_set(const dir_type<Ndim>& dir, typename Set::type data) {
	assert(neighbors[dir] != hpx::invalid_id);
	exchange_semaphore.wait();
	boost::lock_guard<decltype(set_lock)> scope(set_lock);
	(static_cast<Derived*>(this)->*(Set::value))(dir, data);
	exchange_promises[dir].set_value();
	//	set_lock.unlock();
}

template<typename Derived, int Ndim>
template<typename Get, typename Set>
inline void node<Derived, Ndim>::exchange_get(int this_subcycle) {
	using T = typename Get::type;
	wait_my_turn(this_subcycle);
	if (!is_leaf) {
		for (int i = 0; i < Nchild; i++) {
			hpx::apply<action_exchange_get<Get, Set>>(children[i], this_subcycle);
		}
	}
	exchange_promises.clear();
	exchange_promises.resize(Nneighbor);
	std::vector<hpx::shared_future<void>> futures(2 * Nneighbor + 1);
	for (int i = 0; i < Nneighbor; i++) {
		futures[i + Nneighbor] = exchange_promises[i].get_future().share();
	}
	int nn = 0;
	for (dir_type < Ndim > dir; !dir.is_end(); dir++) {
		if (neighbors[dir] != hpx::invalid_id && !dir.is_zero()) {
			auto f = last_operation_future.then(hpx::util::unwrapped([this,dir]() {
				auto tmp = (static_cast<Derived*>(this)->*(Get::value))(dir);
				auto dir_set = dir;
				dir_set.flip();
				return hpx::async<action_exchange_set<Set>>(neighbors[dir], dir_set, tmp);
			}));
			futures[dir] = f.share();
			nn++;
		} else {
			futures[dir] = last_operation_future;
		}
	}
	//	printf("%i %i\n", nn, this->get_self().get_level());
	auto f = last_operation_future.then(hpx::util::unwrapped([this]() {
		for (dir_type<Ndim> dir; !dir.is_end(); dir++) {
			if (neighbors[dir] != hpx::invalid_id && !dir.is_zero()) {
				exchange_semaphore.signal();
			} else {
				exchange_promises[dir].set_value();
			}
		}
	}));
	futures[2 * Nneighbor] = f.share();
	last_operation_future = when_all(futures).share();
	subcycle_promises[subcycle].set_value();
	subcycle++;
}

template<typename Derived, int Ndim>
inline void node<Derived, Ndim>::neighbors_exchange_set(dir_type<Ndim> dir, const std::vector<hpx::id_type>& nephews) {
	assert(neighbors[dir] != hpx::invalid_id);
	exchange_semaphore.wait();
	boost::lock_guard<decltype(set_lock)> scope_lock(set_lock);
	nieces[dir] = nephews;
	exchange_promises[dir].set_value();

}

template<typename Derived, int Ndim>
inline void node<Derived, Ndim>::neighbors_ascend(std::vector<hpx::id_type> input_data) {
	neighbors = std::move(input_data);
	if (!is_leaf) {
		std::vector < std::vector < hpx::id_type >> child_neighbors;
		child_neighbors.resize(Nchild);
		for (child_index_type < Ndim > ci; !ci.is_end(); ci++) {
			child_neighbors[ci].resize(Nneighbor);
			for (dir_type < Ndim > dir; !dir.is_end(); dir++) {
				child_neighbors[ci][dir] = get_sibling_of_child(ci, dir);
			}
		}
		std::vector<hpx::future<void>> cfuts(Nchild);
		for (int i = 0; i < Nchild; i++) {
			cfuts[i] = hpx::async < action_neighbors_ascend > (children[i], child_neighbors[i]);
		}
		wait_all(std::move(cfuts));
	}
}

template<typename Derived, int Ndim>
inline void node<Derived, Ndim>::neighbors_exchange_get(int this_subcycle) {
	wait_my_turn(this_subcycle);
	std::vector<hpx::future<void>> cfuts(Nchild);
	if (!is_leaf) {
		for (int i = 0; i < Nchild; i++) {
			cfuts[i] = hpx::async < node<Derived, Ndim>::action_neighbors_exchange_get > (children[i], this_subcycle);
		}
	}
	exchange_promises.clear();
	exchange_promises.resize(Nneighbor);
	std::vector<hpx::shared_future<void>> futures(2 * Nneighbor + 1);
	for (int i = 0; i < Nneighbor; i++) {
		//	printf("--2 %i\n", this->get_self().get_level());
		futures[i + Nneighbor] = exchange_promises[i].get_future().share();
	}
	int nn = 0;

	dir_type < Ndim > dir;
	dir.set_zero();
	nieces[dir] = children;
	for (dir_type < Ndim > dir; !dir.is_end(); dir++) {
		if (neighbors[dir] != hpx::invalid_id && !dir.is_zero()) {
			auto f =
					last_operation_future.then(
							hpx::util::unwrapped(
									[this,dir]() {
										auto dir_set = dir;
										dir_set.flip();
										return hpx::async<node<Derived, Ndim>::action_neighbors_exchange_set>(neighbors[dir], dir_set, children);
									}));
			futures[dir] = f.share();
			nn++;
		} else {
			futures[dir] = last_operation_future;
		}
	}
	auto f = last_operation_future.then(hpx::util::unwrapped([this]() {
		for (dir_type<Ndim> dir; !dir.is_end(); dir++) {
			if (neighbors[dir] != hpx::invalid_id && !dir.is_zero()) {
				exchange_semaphore.signal();
			} else {
				exchange_promises[dir].set_value();
			}
		}
	}));
	futures[2 * Nneighbor] = f.share();
	last_operation_future = when_all(futures).share();
	subcycle_promises[subcycle].set_value();
	subcycle++;
	if (!is_leaf) {
		wait_all(std::move(cfuts));
	}

}

}

#endif /* NODE_HPP_ */
