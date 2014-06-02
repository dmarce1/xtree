/*
 * tree.hpp
 *
 *  Created on: May 24, 2014
 *      Author: dmarce1
 */

#ifndef TREE_HPP_
#define TREE_HPP_

namespace xtree {

template<typename Member, int Ndim>
class tree: public hpx::components::managed_component_base<tree<Member, Ndim>, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr> {
public:
	static constexpr int Nbranch = 2;
	static constexpr int Nneighbor = pow_<3, Ndim>::value;
	const char* name = "tree";
private:
	using base_type = hpx::components::managed_component_base<tree<Member,Ndim>, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr>;
	using component_type = hpx::components::managed_component<tree<Member,Ndim>>;
	using silo_output_type = silo_output<Ndim, 1>;
	hpx::id_type silo_gid;
	hpx::id_type this_gid;
	hpx::id_type root_node_gid;
	std::set<node<Member, Ndim>*> nodes;
	tree<Member, Ndim>* this_ptr;
	std::array<hpx::id_type, Nbranch> child_gids;
	load_balancer* load_balancer_ptr;
	hpx::id_type load_balancer_gid;
	hpx::lcos::local::mutex dir_lock;
public:

	tree() {
		assert(false);
	}
	tree(component_type* back_ptr) :
			base_type(back_ptr) {
		static bool initialized = false;
		assert(!initialized);
		initialized = true;
		int my_id, id_cnt;
		auto localities = hpx::find_all_localities();
		std::vector < hpx::future < hpx::id_type >> futures(Nbranch);

		my_id = hpx::get_locality_id();
		id_cnt = localities.size();
		this_gid = base_type::get_gid();
		this_ptr = this;
		hpx::register_id_with_basename(name, this_gid, my_id);

		for (int i = 0; i < Nbranch; i++) {
			int j = my_id * Nbranch + i + 1;
			if (j < localities.size()) {
				futures[i] = hpx::new_<tree<Member, Ndim>>(localities[j]);
			} else {
				futures[i] = hpx::make_ready_future(hpx::invalid_id);
			}
		}
		for (int i = 0; i < Nbranch; i++) {
			child_gids[i] = futures[i].get();
		}
		load_balancer_gid = hpx::new_ < load_balancer > (hpx::find_here()).get();
		load_balancer_ptr = (hpx::async<load_balancer::action_get_ptr>(load_balancer_gid)).get();
		if (my_id == 0) {
			hpx::future < hpx::id_type > fut1;
			hpx::future < hpx::id_type > fut2;
			fut1 = hpx::new_ < silo_output_type > (hpx::find_here());
			fut2 = new_node(location<Ndim>(), hpx::invalid_id, vector<hpx::id_type, Nneighbor>(hpx::invalid_id));
			silo_gid = fut1.get();
			root_node_gid = fut2.get();
		} else {
			silo_gid = hpx::invalid_id;
			root_node_gid = hpx::invalid_id;
		}
	}
	template<typename Archive>
	void serialize(Archive& ar, const int v) {
		assert(false);
	}

	virtual ~tree() {
	}

	hpx::future<hpx::id_type> get_new_node(const location<Ndim>& _loc, hpx::id_type _parent_id, const vector<hpx::id_type, Nneighbor>& _neighbors) {
		hpx::shared_future < hpx::id_type > id_future = hpx::new_<node<Member, Ndim>>(hpx::find_here(), _loc, _parent_id, std::move(_neighbors), this).share();
		auto fut1 = id_future.then(hpx::util::unwrapped([](hpx::id_type id) {
			return hpx::async<typename node<Member,Ndim>::action_get_this>(id);
		}));
		return fut1.then(hpx::util::unwrapped([this,id_future](node<Member,Ndim>* ptr) {
			dir_lock.lock();
			auto test = nodes.insert(ptr);
			dir_lock.unlock();
			assert(test.second);
			return id_future;
		}));
	}

	hpx::future<hpx::id_type> new_node(const location<Ndim>& _loc, hpx::id_type _parent_id, const vector<hpx::id_type, Nneighbor>& _neighbors) {
		return load_balancer_ptr->increment_load().then(hpx::util::unwrapped([_loc,_parent_id,_neighbors](hpx::id_type id) {
			return hpx::async<action_get_new_node>(id, _loc,_parent_id,_neighbors);
		}));
	}

	void delete_node(node<Member, Ndim>* ptr) {
		load_balancer_ptr->decrement_load();
		dir_lock.lock();
		nodes.erase(nodes.find(ptr));
		dir_lock.unlock();
	}

	template<typename Op>
	void begin_execute() {
		hpx::future<bool> fut;
		int lev;
		switch (Op::op) {
		case op_type::AMR_ASCEND:
			printf("Executing AMR_ASCEND\n");
			lev = 0;
			while (this_ptr->execute<Op>(lev).get()) {
				lev++;
			}
			break;
		case op_type::ASCEND:
			printf("Executing ascend\n");
			lev = 0;
			while (this_ptr->execute<Op>(lev).get()) {
				lev++;
			}
			break;
		case op_type::DESCEND:
			printf("Executing descend\n");
			lev = get_max_level().get();
			while (lev >= 0) {
				this_ptr->execute<Op>(lev).get();
				lev--;
			}
			break;
		case op_type::EXCHANGE:
		case op_type::REBRANCH:
			this_ptr->execute<Op>(-1).get();
			break;
		default:
			assert(false);
			break;
		}

	}

	hpx::future<int> get_max_level() {
		int maxlev = 0;
		std::vector<hpx::future<int>> futures(Nbranch);
		for (int i = 0; i < Nbranch; i++) {
			if (child_gids[i] != hpx::invalid_id) {
				futures[i] = hpx::async<action_get_max_level>(child_gids[i]);
			} else {
				futures[i] = hpx::make_ready_future(0);
			}
		}
		for (auto i = nodes.begin(); i != nodes.end(); i++) {
			maxlev = std::max((*i)->get_level(), maxlev);
		}
		return when_all(futures).then(hpx::util::unwrapped([maxlev](std::vector<hpx::future<int>> futures) {

			int _maxlev = maxlev;
			for( int i = 0; i < Nbranch; i++) {
				_maxlev = std::max(_maxlev,futures[i].get());
			}
			return _maxlev;

		}));
	}

	hpx::shared_future<void> get_terminal_future() {
		std::vector<hpx::shared_future<void>> futures(Nbranch + nodes.size());
		for (int i = 0; i < Nbranch; i++) {
			if (child_gids[i] != hpx::invalid_id) {
				futures[i] = hpx::async<action_get_terminal_future>(child_gids[i]);
			} else {
				futures[i] = hpx::make_ready_future();
			}
		}
		int j = 0;
		for (auto i = nodes.begin(); i != nodes.end(); i++) {
			futures[Nbranch + j++] = (*i)->get_last_future();
		}
		return when_all(futures).share();
	}
	template<typename Ops, int Iter = std::tuple_size<Ops>::value>
	struct execute_ops {
		static hpx::shared_future<void> run(tree* _this) {
			constexpr int ThisIter = std::tuple_size < Ops > ::value - Iter;
			printf("Executing op %i\n", ThisIter);
			_this->begin_execute<typename std::tuple_element<ThisIter, Ops>::type>();
			return execute_ops<Ops, Iter - 1>::run(_this);
		}
	};
	template<typename Ops>
	struct execute_ops<Ops, 1> {
		static hpx::shared_future<void> run(tree* _this) {
			constexpr int ThisIter = std::tuple_size < Ops > ::value - 1;
			printf("Executing op %i\n", ThisIter);
			_this->begin_execute<typename std::tuple_element<ThisIter, Ops>::type>();
			return _this->get_terminal_future();
		}
	};
	template<typename Ops>
	hpx::shared_future<void> execute_operators() {
		dir_lock.lock();
		printf("Executing operators\n");
		auto future = execute_ops<Ops, std::tuple_size<Ops>::value>::run(this);
		dir_lock.unlock();
		return future;
	}

	template<typename Op>
	hpx::future<bool> execute(int level) {
		bool rc = false;
		using action = action_execute<Op>;
		std::vector<hpx::future<bool>> futures(Nbranch);
		for (int i = 0; i < Nbranch; i++) {
			if (child_gids[i] != hpx::invalid_id) {
				futures[i] = hpx::async<action>(child_gids[i], level);
			} else {
				futures[i] = hpx::make_ready_future(false);
			}
		}
		for (auto i = nodes.begin(); i != nodes.end(); i++) {
			if (level == -1 || ((*i)->get_level() == level)) {
				rc = true;
				(*i)->template setup_op_dataflow<Op>();
			}
		}
		return when_all(futures).then(hpx::util::unwrapped([rc](std::vector<hpx::future<bool>> futures) {
			bool _rc = false;
			if( !rc ) {
				for( int i = 0; i < Nbranch; i++) {
					if( futures[i].get() ) {
						_rc = true;
						break;
					}
				}
			}
			return rc || _rc;
		}));
	}

	template<op_type Op, typename T, typename Get, typename Set>
	struct operation {
		static constexpr op_type op = Op;
		const Get get;
		const Set set;
		typedef T type;
		operation() :
				get(), set() {
		}
	};

	template<typename Op>
	using action_execute = typename hpx::actions::make_action<decltype(&tree::execute<Op>), &tree::execute<Op>>::type;
	template<typename Ops>
	using action_execute_operators = typename hpx::actions::make_action<decltype(&tree::execute_operators<Ops>), &tree::execute_operators<Ops>>::type;XTREE_MAKE_ACTION( action_get_new_node, tree::get_new_node );XTREE_MAKE_ACTION( action_get_max_level, tree::get_max_level );XTREE_MAKE_ACTION( action_get_terminal_future, tree::get_terminal_future );

};

}

#endif /* TREE_HPP_ */
