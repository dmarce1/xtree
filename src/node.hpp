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

template<typename Derived, int Ndim>
class node: public hpx::components::abstract_managed_component_base<
		node<Derived, Ndim>> {
public:
	static constexpr int Nchild = pow_<2, Ndim>::value;
	static constexpr int Nneighbor = pow_<3, Ndim>::value;
	using base_type = hpx::components::managed_component_base<Derived>;
	using child_array_type = std::vector<hpx::id_type>;
	using niece_array_type = std::vector<std::vector<hpx::id_type>>;
	using neighbor_array_type = std::vector<hpx::id_type>;
	using wrapped_type = Derived;

private:
	static bool child_is_niece_of(child_index_type<Ndim> ci,
			dir_type<Ndim> dir) {
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
//	hpx::lcos::local::mutex subcycle_lock;
	mutable hpx::lcos::local::mutex branch_lock;
	mutable hpx::lcos::local::mutex sem_lock;
	mutable hpx::lcos::local::mutex set_lock;
private:
	hpx::id_type get_sibling_of_child(child_index_type<Ndim> ci,
			dir_type<Ndim> dir) {
		child_index_type<Ndim> nci;
		dir_type<Ndim> ndi;
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
public:
	node() :
			exchange_semaphore(0) {
		//	printf( "loading node\n");
	}
	wrapped_type* initialize(const location<Ndim>& _loc,
			hpx::util::tuple<hpx::id_type, neighbor_array_type> other_ids,
			tree<wrapped_type, Ndim>* ltree) {
		subcycle = 0;
		last_operation_future = hpx::make_ready_future();
		local_tree = ltree;
		self = _loc;
		//	cur_fut = prev_fut = hpx::make_ready_future();
		neighbors = hpx::util::get < 1 > (other_ids);
		is_leaf = true;
		is_branching = false;
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
		parent = hpx::util::get < 0 > (other_ids);
		static_cast<Derived*>(this)->initialize();
		return dynamic_cast<wrapped_type*>(this);
	}

	virtual ~node() {
		if (!is_leaf) {
			debranch().get();
		}
		if (get_level() != 0) {
			local_tree->delete_node(static_cast<wrapped_type*>(this));
		}
		parent = hpx::invalid_id;
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

	hpx::shared_future<void> branch() {
		if (!is_leaf) {
			return hpx::make_ready_future();
		}
		is_branching = true;
		std::vector < hpx::future < hpx::id_type >> futs(Nchild);
		for (child_index_type<Ndim> ci; !ci.is_end(); ci++) {
			neighbor_array_type pack(Nneighbor);
			//	for (dir_type<Ndim> dir; !dir.is_end(); dir++) {
			//	pack[dir] = get_sibling_of_child(ci, dir);
			//	}
			auto tgid = static_cast<Derived*>(this)->get_gid();
			//		printf("Making new child\n");
			futs[ci] = local_tree->new_node(self.get_child(ci), tgid, pack, 0);
		}
		hpx::future<void> f =
				when_all(futs).then(
						hpx::util::unwrapped(
								[this](std::vector < hpx::future < hpx::id_type >> futs) {
									is_leaf = false;
									dir_type<Ndim> zero;
									zero.set_zero();
									for (int ci = 0; ci < Nchild; ci++) {
										children[ci] = futs[ci].get();
										nieces[zero][ci] = children[ci];
									}
								}));
		return f.share();
	}

	hpx::future<void> debranch() {
		branch_lock.lock();
		is_leaf = true;
		std::vector<hpx::future<void>> nfutures(Nneighbor);
		std::vector<hpx::future<void>> cfutures(Nchild);
		for (child_index_type<Ndim> ci; !ci.is_end(); ci++) {
			if (children[ci] != hpx::invalid_id) {
				cfutures[ci] = hpx::async < action_debranch > (children[ci]);
			} else {
				cfutures[ci] = hpx::make_ready_future();
			}
		}
		/*		for (dir_type<Ndim> dir; !dir.is_end(); dir++) {
		 if (neighbors[dir] != hpx::invalid_id) {
		 nfutures[dir] = hpx::async<action_notify_debranch>(neighbors[dir], dir);
		 } else {
		 nfutures[dir] = hpx::make_ready_future();
		 }
		 }*/
		auto fut1 = when_all(cfutures).then(
				hpx::util::unwrapped([this](std::vector<hpx::future<void>>) {
					branch_lock.lock();
					for (std::size_t ci = 0; ci != Nchild; ++ci) {
						children[ci] = hpx::invalid_id;
					}
					branch_lock.unlock();
				}));
		branch_lock.unlock();
		return fut1;
	}
	int get_level() const {
		return self.get_level();
	}

	int get_subcycle() const {
		return subcycle;
	}
	wrapped_type* get_this() {
		return static_cast<wrapped_type*>(this);
	}

	bound_type get_boundary_type(const dir_type<Ndim>& dir) const {
		if (neighbors[dir] != hpx::invalid_id) {
			return DECOMP;
		} else {
			for (int di = 0; di < Ndim; di++) {
				if (((dir[di] == -1) && (self.get_location(di) == 0))
						|| ((dir[di] == +1)
								&& (self.get_location(di)
										== ((1 << self.get_level()) - 1)))) {
					return PHYS;
				}
			}
			return AMR;
		}
	}

	bool has_amr_boundary() const {
		for (dir_type<Ndim> dir; !dir.is_end(); dir++) {
			if (get_boundary_type(dir) == AMR) {
				return true;
			}
		}
		return false;
	}

	void wait_my_turn(int this_subcycle) {
		//	subcycle_lock.lock();
		while (this_subcycle != subcycle) {
			//	subcycle_lock.unlock();
			hpx::this_thread::suspend();
			//		subcycle_lock.lock();
		}
		//	subcycle_lock.unlock();
	}

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
		virtual ~operation_base() {
		}
		virtual void operator()(node&) const = 0;
		virtual bool is_regrid() const {
			return false;
		}
	};

	using operation_type = std::shared_ptr<operation_base>;

	hpx::shared_future<void> operations_end(int this_subcycle) {
		wait_my_turn(this_subcycle);
		hpx::future<void> future;
		if (!is_leaf) {
			std::vector<hpx::future<void>> futs(Nchild);
			for (std::size_t i = 0; i != Nchild; ++i) {
				futs[i] = hpx::async < action_operations_end
						> (children[i], this_subcycle);
			}
			future = when_all(futs);
		} else {
			future = hpx::make_ready_future();
		}
		auto f =
				when_all(future, last_operation_future).then(
						hpx::util::unwrapped(
								[this](hpx::util::tuple<hpx::future<void>, hpx::shared_future<void>>) {
									subcycle = 0;
								}));
		last_operation_future = f.share();
		printf("%i\n", this->get_self().get_level());

		return last_operation_future;
	}

	void execute_operations(std::vector<operation_type> operations) {
		hpx::shared_future<void> fut;
		int i;
		fut = hpx::make_ready_future().share();
		for (i = 0; i < operations.size(); i++) {
			auto f = fut.then(hpx::util::unwrapped([i,this,&operations]() {
				operations[i]->operator()(*this);
			}));
			fut = f.share();
		}
		//		printf("Done with ops - pre----1\n");
		fut.get();
		if (!operations[i - 1]->is_regrid()) {
			fut = operations_end(subcycle);
		} else {
			fut = last_operation_future;
		}
		//		printf("Done with ops - pre----2\n");
		//		printf("Done with ops - pre\n");
		fut.get();
		//		printf("Done with ops\n");
	}

	template<typename Function>
	struct local_operation: operation_base {
		void operator()(node& root) const {
			root.local<Function>(root.get_subcycle());
		}
	};

	template<typename Function>
	struct regrid_operation: operation_base {
		virtual bool is_regrid() const {
			return true;
		}
		void operator()(node& root) const {
			root.regrid<Function>(root.get_subcycle());
		}
	};

	template<typename Function>
	struct ascend_operation: operation_base {
		void operator()(node& root) const {
			root.ascend<Function>(
					hpx::make_ready_future<typename Function::type>(
							typename Function::type()), root.get_subcycle());
		}
	};

	template<typename Function>
	struct descend_operation: operation_base {
		void operator()(node& root) const {
			root.descend<Function>(root.get_subcycle());
		}
	};

	template<typename Get, typename Set>
	struct exchange_operation: operation_base {
		void operator()(node& root) const {
			root.exchange_get<Get, Set>(root.get_subcycle());
		}
	};

	template<typename Function>
	static operation_type make_local_operation() {
		auto operation = std::make_shared<local_operation<Function>>();
		return std::static_pointer_cast < operation_base > (operation);
	}

	template<typename Function>
	static operation_type make_regrid_operation() {
		auto operation = std::make_shared<regrid_operation<Function>>();
		return std::static_pointer_cast < operation_base > (operation);
	}

	template<typename Function>
	static operation_type make_ascend_operation() {
		auto operation = std::make_shared<ascend_operation<Function>>();
		return std::static_pointer_cast < operation_base > (operation);
	}

	template<typename Function>
	static operation_type make_descend_operation() {
		auto operation = std::make_shared<descend_operation<Function>>();
		return std::static_pointer_cast < operation_base > (operation);
	}

	template<typename Get, typename Set>
	static operation_type make_exchange_operation() {
		auto operation = std::make_shared<exchange_operation<Get, Set>>();
		return std::static_pointer_cast < operation_base > (operation);
	}

	template<typename Function>
	void local(int this_subcycle) {
		wait_my_turn(this_subcycle);
		std::vector<hpx::future<void>> futures;
		hpx::shared_future<void> rc;
		if (!is_leaf) {
			futures.resize(Nchild);
			for (int i = 0; i < Nchild; ++i) {
				futures[i] = hpx::async < action_local
						< Function >> (children[i], this_subcycle);
			}
		}
		rc = hpx::async([this]() {
			(static_cast<Derived*>(this)->*(Function::value))();
		}).share();

		if (futures.size()) {
			rc = when_all(rc, when_all(futures)).share();
		}

		//	subcycle_lock.lock();
		last_operation_future = rc;
		subcycle++;
//		subcycle_lock.unlock();
	}

	template<typename Function>
	bool regrid(int this_subcycle) {
		//		printf("Regridding\n");
		//		wait_my_turn(this_subcycle);
		std::vector<hpx::future<bool>> futures;
		bool rc;
		if (!is_leaf) {
			futures.resize(Nchild);
			for (int i = 0; i < Nchild; ++i) {
				auto c = children[i];
				futures[i] = hpx::async < action_regrid
						< Function >> (c, this_subcycle);
			}
			wait_all(futures);
			for (auto i = 0; i != Nchild; i++) {
				if (!futures[i].get()) {
					debranch().get();
				}
			}
			printf("Inside\n");
			rc = true;
		} else {
			rc = (static_cast<Derived*>(this)->*(Function::value))();
			if (rc) {
				branch().get();
				printf("Ready\n");
			}
		}

		if (this->get_self().get_level() == 0) {

//			hpx::async < action_find_neighbors
//					> (static_cast<Derived*>(this)->get_gid()).get();
			//		subcycle = 0;
			neighbors_exchange_get(subcycle);
			neighbors_ascend(
					hpx::make_ready_future < std::vector<hpx::id_type>
							> (std::vector < hpx::id_type > (Nneighbor)),
					subcycle);
			operations_end(subcycle).get();
		}

		//	subcycle_lock.lock();
//		last_operation_future = hpx::make_ready_future();
		//	subcycle++;
		//	subcycle_lock.unlock();
		//		printf("Done Regridding\n");
		return rc;
	}

	template<typename Function>
	void ascend(hpx::future<typename Function::type> input_data_future,
			int this_subcycle) {
		using T = typename Function::type;
		wait_my_turn(this_subcycle);
		auto promises = std::make_shared<std::vector<hpx::promise<T>>>();
		promises->resize(Nchild);
		auto f =
				when_all(input_data_future, last_operation_future).then(
						hpx::util::unwrapped(
								[=](HPX_STD_TUPLE<hpx::future<T>,hpx::shared_future<void>> tupfut) {
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
				hpx::apply < action_ascend
						< Function
								>> (children[i], (*promises)[i].get_future(), this_subcycle);
			}
		}
		//	subcycle_lock.lock();
		last_operation_future = f.share();
		subcycle++;
		//	subcycle_lock.unlock();
	}

	template<typename Function>
	hpx::shared_future<typename Function::type> descend(int this_subcycle) {
		//test = 0;
		using T = typename Function::type;
		wait_my_turn(this_subcycle);
		hpx::shared_future<T> return_future;
		hpx::future<T> f;
		if (!is_leaf) {
			std::vector < hpx::future < T >> futures(Nchild);
			for (int i = 0; i < Nchild; i++) {
				futures[i] = hpx::async < action_descend
						< Function >> (children[i], this_subcycle);
			}
			auto f =
					when_all(when_all(futures), last_operation_future).then(
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

			auto f =
					last_operation_future.then(
							hpx::util::unwrapped(
									[this]() {
										std::vector<T> empty_vector;
										return (static_cast<Derived*>(this)->*(Function::value))(empty_vector);
									}));

			return_future = f.share();
		}
//		subcycle_lock.lock();
		last_operation_future = return_future;
		subcycle++;
//		subcycle_lock.unlock();
		return return_future;
	}

	template<typename Set>
	void exchange_set(const dir_type<Ndim>& dir, typename Set::type data) {
		set_lock.lock();
		assert(neighbors[dir] != hpx::invalid_id);
		exchange_semaphore.wait();
		(static_cast<Derived*>(this)->*(Set::value))(dir, data);
		exchange_promises[dir].set_value();
		set_lock.unlock();
	}

	template<typename Get, typename Set>
	void exchange_get(int this_subcycle) {
		using T = typename Get::type;
		wait_my_turn(this_subcycle);
		if (!is_leaf) {
			for (int i = 0; i < Nchild; i++) {
				hpx::apply<action_exchange_get<Get, Set>>(children[i],
						this_subcycle);
			}
		}
		//	subcycle_lock.lock();
		exchange_promises.clear();
		exchange_promises.resize(Nneighbor);
		std::vector<hpx::shared_future<void>> futures(2 * Nneighbor + 1);
		for (int i = 0; i < Nneighbor; i++) {
			futures[i + Nneighbor] = exchange_promises[i].get_future();
		}
		int nn = 0;
		for (dir_type<Ndim> dir; !dir.is_end(); dir++) {
			if (neighbors[dir] != hpx::invalid_id && !dir.is_zero()) {
				auto f =
						last_operation_future.then(
								hpx::util::unwrapped(
										[this,dir]() {
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
		subcycle++;
//		subcycle_lock.unlock();

	}

	template<typename Function>
	using action_local = hpx::actions::make_action<void(node::*)(int), &node::local<Function>>;
//
	template<typename Function>
	using action_regrid = hpx::actions::make_action< bool(node::*)(int), &node::regrid<Function>>;
//
	template<typename Function>
	using action_ascend = hpx::actions::make_action< void(node::*)(hpx::future<typename Function::type>, int), &node::ascend<Function>>;
//
	template<typename Function>
	using action_descend = hpx::actions::make_action<hpx::shared_future<typename Function::type>(node::*)(int), &node::descend<Function>>;
//
	template<typename Set>
	using action_exchange_set = hpx::actions::make_action<void(node::*)(const dir_type<Ndim>&, typename Set::type), &node::exchange_set<Set>>;
//
	template<typename Get, typename Set>
	using action_exchange_get = hpx::actions::make_action<void(node::*)(int), &node::exchange_get<Get,Set>>;

	void neighbors_exchange_set(dir_type<Ndim> dir,
			const std::vector<hpx::id_type>& nephews) {
		set_lock.lock();
		assert(neighbors[dir] != hpx::invalid_id);
		exchange_semaphore.wait();
		nieces[dir] = nephews;
		exchange_promises[dir].set_value();
		set_lock.unlock();

	}

	void neighbors_ascend(
			hpx::future<std::vector<hpx::id_type>> input_data_future,
			int this_subcycle);
	void neighbors_exchange_get(int this_subcycle);
	using action_neighbors_ascend = hpx::actions::make_action<void(node::*)(hpx::future<std::vector<hpx::id_type>>,
			int), &node::neighbors_ascend>;
	using action_neighbors_exchange_get = hpx::actions::make_action<void(node::*)(int), &node::neighbors_exchange_get>;
	using action_neighbors_exchange_set = hpx::actions::make_action<void(node::*)(dir_type<Ndim>,
			const std::vector<hpx::id_type>&), &node::neighbors_exchange_set>;
//
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, operations_end, action_operations_end);//
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, initialize, action_initialize);	//
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, debranch, action_debranch);	//
	HPX_DEFINE_COMPONENT_ACTION_TPL(node, get_this, action_get_this);	//
//
//
}
;

template<typename Derived, int Ndim>
void node<Derived, Ndim>::neighbors_ascend(
		hpx::future<std::vector<hpx::id_type>> input_data_future,
		int this_subcycle) {
	wait_my_turn(this_subcycle);
	auto promises = std::make_shared<
			std::vector<hpx::promise<std::vector<hpx::id_type>>> >();
	promises->resize(Nchild);
	auto f =
			when_all(input_data_future, last_operation_future).then(
					hpx::util::unwrapped(
							[=](HPX_STD_TUPLE<hpx::future<std::vector<hpx::id_type>>,hpx::shared_future<void>> tupfut) {
								auto& g = HPX_STD_GET(0, tupfut);
								neighbors = g.get();
								if( !is_leaf ) {
									std::vector < std::vector < hpx::id_type >> child_neighbors;
									child_neighbors.resize(Nchild);
									for (child_index_type<Ndim> ci; !ci.is_end(); ci++) {
										child_neighbors[ci].resize(Nneighbor);
										for (dir_type<Ndim> dir; !dir.is_end(); dir++) {
											if( !dir.is_zero()) {
												child_neighbors[ci][dir] = get_sibling_of_child(ci, dir);
											}
										}
									}
									for( int i = 0; i < Nchild; i++) {
										(*promises)[i].set_value(child_neighbors[i]);
									}
								}
							}));
	if (!is_leaf) {
		for (int i = 0; i < Nchild; i++) {
			printf("--1\n");
			hpx::apply < action_neighbors_ascend
					> (children[i], (*promises)[i].get_future(), this_subcycle);
		}
	}
//	subcycle_lock.lock();
	last_operation_future = f.share();
	subcycle++;
//	subcycle_lock.unlock();
}

template<typename Derived, int Ndim>
void node<Derived, Ndim>::neighbors_exchange_get(int this_subcycle) {
	wait_my_turn(this_subcycle);
	if (!is_leaf) {
		for (int i = 0; i < Nchild; i++) {
			hpx::apply < node<Derived, Ndim>::action_neighbors_exchange_get
					> (children[i], this_subcycle);
		}
	}
//	subcycle_lock.lock();
	exchange_promises.clear();
	exchange_promises.resize(Nneighbor);
	std::vector<hpx::shared_future<void>> futures(2 * Nneighbor + 1);
	for (int i = 0; i < Nneighbor; i++) {
		printf("--2 %i\n", this->get_self().get_level());
		futures[i + Nneighbor] = exchange_promises[i].get_future();
	}
	int nn = 0;

	dir_type<Ndim> dir;
	dir.set_zero();
	nieces[dir] = children;
	for (dir_type<Ndim> dir; !dir.is_end(); dir++) {
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
	subcycle++;

}

}

#endif /* NODE_HPP_ */
