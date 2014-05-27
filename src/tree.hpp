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
class tree : public hpx::components::managed_component_base<tree<Member,Ndim>>  {
public:

	static constexpr int Nbranch = 2;
	static constexpr int Nneighbor = pow_<3, Ndim>::get();

private:
	typedef hpx::components::managed_component_base<tree<Member,Ndim>> base_type;
	static const std::pair<int, hpx::id_type> mins_begin;
	static std::vector<hpx::id_type> neighbors;
	static hpx::lcos::local::counting_semaphore semaphore;
	static int my_load;
	static std::set<node<Member,Ndim>*> nodes;
	static vector<hpx::id_type, Nbranch> child_ids;
	static vector<hpx::id_type, Nbranch> child_gids;
	static hpx::id_type this_gid;
	static tree<Member,Ndim>* this_ptr;

public:

	tree() {
		this_gid = base_type::get_gid();
		this_ptr = this;
		std::vector<hpx::future<hpx::id_type>> futures(Nbranch);
		auto localities = hpx::find_all_localities();
		int mask, neighbor_id, my_id, id_cnt;

		my_id = hpx::get_locality_id();
		id_cnt = localities.size();

		for (int i = 0; i < Nbranch; i++) {
			int j = my_id * Nbranch + i + 1;
			if (j < localities.size()) {
				child_ids[i] = localities[j];
				futures[i] = hpx::new_<tree<Member,Ndim>>(localities[j]);
			} else {
				child_ids[i]  = hpx::invalid_id;
				futures[i] = hpx::make_ready_future(hpx::invalid_id);
			}
		}

		std::list<int> tmp;
		tmp.push_front(my_id);
		for (mask = 1; mask > 0; mask <<= 1) {
			neighbor_id = mask ^ my_id;
			if (neighbor_id < id_cnt) {
				tmp.push_front(neighbor_id);
			}
		}
		std::vector<int> tmp2 = std::vector<int>(tmp.begin(), tmp.end());
		std::sort(tmp2.begin(), tmp2.end());
		neighbors.resize(tmp2.size());
		for (size_t i = 0; i < neighbors.size(); i++) {
			neighbors[i] =localities[tmp2[i]];
		}

		hpx::lcos::wait_all(futures);
		for( int i = 0; i < Nbranch; i++) {
			child_gids[i] = futures[i].get();
		}
		if (my_id == 0) {
			new_node(location<Ndim>(), hpx::invalid_id, vector<hpx::id_type, Nneighbor>(hpx::invalid_id));
		}

	}

	virtual ~tree() {
	}

	static hpx::future<hpx::id_type>  get_new_node(const location<Ndim>& _loc, hpx::id_type _parent_id, vector<hpx::id_type, Nneighbor> _neighbors) {
		hpx::shared_future<hpx::id_type> id_future = (hpx::new_<node<Member,Ndim>>(hpx::find_here(), _loc,_parent_id,_neighbors)).share();
		auto fut1 = id_future.then(hpx::util::unwrapped([](hpx::id_type id) {
			return hpx::async<typename node<Member,Ndim>::action_get_this>(id);
		}));
		return fut1.then(hpx::util::unwrapped([id_future](node<Member,Ndim>* ptr){
			nodes.insert(ptr);
			return id_future;
		}));
	}


	static hpx::future<hpx::id_type> new_node(const location<Ndim>& _loc, hpx::id_type _parent_id, vector<hpx::id_type, Nneighbor> _neighbors) {
		return increment_load().then(hpx::util::unwrapped([_loc,_parent_id,_neighbors](hpx::id_type id){
			return hpx::async<action_get_new_node>(id, _loc,_parent_id,_neighbors);
		}));
	}

	static void delete_node(node<Member,Ndim>* ptr) {
		nodes.erase(nodes.find(ptr));
	}

	template<typename Op>
	static void begin_execute() {
		hpx::future<bool> fut;
		int lev;
		switch (Op::op) {
		case op_type::ASCEND:
			lev = 0;
			while (this_ptr->execute<Op>(lev).get()) {
				lev++;
			}
			break;
		case op_type::DESCEND:
			lev = get_max_level().get();
			while (lev >= 0) {
				this_ptr->execute<Op>(lev).get();
				lev--;
			}
			break;
		case op_type::EXCHANGE:
			this_ptr->execute<Op>(-1).get();
			break;
		}

	}

	template<typename Op>
	hpx::future<bool> execute(int level) {
		bool rc = false;
		using action = action_execute<Op>;
		std::vector<hpx::future<bool>> futures(Nbranch + nodes.size());
		for (int i = 0; i < Nbranch; i++) {
			if (child_ids[i] != hpx::invalid_id) {
				futures[i] = hpx::async<action>(child_gids[i], level);
			} else {
				futures[i] = hpx::make_ready_future(true);
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

	static hpx::future<int> get_max_level() {
		int maxlev = 0;
		std::vector<hpx::future<int>> futures(Nbranch);
		for (int i = 0; i < Nbranch; i++) {
			if (child_ids[i] != hpx::invalid_id) {
				futures[i] = hpx::async<action_get_max_level>(child_ids[i]);
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


	static hpx::shared_future<void> get_terminal_future() {
		std::vector<hpx::shared_future<void>> futures(Nbranch + nodes.size());
		for (int i = 0; i < Nbranch; i++) {
			if (child_ids[i] != hpx::invalid_id) {
				futures[i] = hpx::async<action_get_terminal_future>(child_ids[i]);
			} else {
				futures[i] = hpx::make_ready_future();
			}
		}
		int j = 0;
		for (auto i = nodes.begin(); i != nodes.end(); i++) {
			futures[Nbranch + j++] =(*i)->get_last_future();
		}
		return when_all(futures).share();
	}

	static hpx::future<std::pair<int, hpx::id_type>> lock_servlet(std::pair<int, hpx::id_type> best_mins, std::list<hpx::id_type> remaining) {
		hpx::future<std::pair<int, hpx::id_type>> fut;

		semaphore.wait();
		if (my_load < best_mins.first) {
			if (best_mins.second != hpx::invalid_id) {
				hpx::apply<action_unlock_servlet>(best_mins.second, false);
			}
			best_mins.first = my_load;
			best_mins.second = hpx::find_here();
		} else {
			semaphore.signal();
		}

		if (remaining.size() > 0) {
			const auto id = remaining.front();
			remaining.pop_front();
			fut = hpx::async<action_lock_servlet>(id, best_mins, remaining);
		} else {
			fut = hpx::make_ready_future(best_mins);
		}

		return fut;
	}

	static int get_load() {
		int i;
		semaphore.wait();
		i = my_load;
		semaphore.signal();
		return i;
	}

	static hpx::future<hpx::id_type> increment_load() {

		std::pair<int, hpx::id_type> mins;
		mins.first = INT_MAX;
		mins.second = hpx::invalid_id;
		std::list<hpx::id_type> remaining;
		for (size_t i = 1; i < neighbors.size(); i++) {
			remaining.push_back(neighbors[i]);
		}
		auto fut1 = hpx::async<action_lock_servlet>(neighbors[0], mins_begin, remaining);
		auto fut2 = fut1.then(hpx::util::unwrapped([](std::pair<int, hpx::id_type> mins) {
			return hpx::async<action_unlock_servlet>(mins.second, true);
		}));
		return fut2;
	}

	static void decrement_load() {
		semaphore.wait();
		my_load--;
		semaphore.signal();
	}

	static hpx::id_type unlock_servlet(bool inc_cnt) {

		if (inc_cnt) {
			my_load++;
		}
		semaphore.signal();
		return hpx::find_here();
	}


	template<op_type Op, typename T, get_type<Member, T, Ndim> Get, set_type<Member, T, Ndim> Set = nullptr>
	struct operation {
		static constexpr op_type op = Op;
		static constexpr get_type<Member, T, Ndim> get = Get;
		static constexpr set_type<Member, T, Ndim> set = Set;
		typedef T type;
	};



private:
	template<typename Ops, int Iter >
	struct execute_ops {
		static hpx::shared_future<void> run(tree* _this) {
			_this->begin_execute<typename std::tuple_element<std::tuple_size<Ops>::value-Iter,Ops>::type>();
			return execute_ops<Ops,Iter-1>::run(_this);
		}
	};

	template<typename Ops>
	struct execute_ops<Ops,1> {
		static hpx::shared_future<void> run(tree* _this) {
			_this->begin_execute<typename std::tuple_element<std::tuple_size<Ops>::value-1,Ops>::type>();
			return _this->get_terminal_future();
		}
	};

public:
	template<typename Ops>
	hpx::shared_future<void> execute_operators() {
		return execute_ops<Ops, std::tuple_size<Ops>::value >::run(this);
	}


	template<typename Op>
	using action_execute = typename hpx::actions::make_action<decltype(&tree::execute<Op>), &tree::execute<Op>>::type;

	XTREE_MAKE_ACTION( action_get_new_node, tree::get_new_node );
	XTREE_MAKE_ACTION( action_get_max_level, tree::get_max_level );
	XTREE_MAKE_ACTION( action_get_terminal_future, tree::get_terminal_future );
	XTREE_MAKE_ACTION( action_lock_servlet, tree::lock_servlet );
	XTREE_MAKE_ACTION( action_unlock_servlet, tree::unlock_servlet );


};

template<typename Member, int Ndim>
const std::pair<int, hpx::id_type> tree<Member,Ndim>::mins_begin = std::pair<int, hpx::id_type>(INT_MAX, hpx::invalid_id);


template<typename Member, int Ndim>
std::vector<hpx::id_type> tree<Member,Ndim>::neighbors;

template<typename Member, int Ndim>
hpx::lcos::local::counting_semaphore tree<Member,Ndim>::semaphore;


template<typename Member, int Ndim>
int tree<Member,Ndim>::my_load;


template<typename Member, int Ndim>
std::set<node<Member,Ndim>*> tree<Member,Ndim>::nodes;

template<typename Member, int Ndim>
vector<hpx::id_type,tree<Member,Ndim>::Nbranch> tree<Member,Ndim>::child_ids;

template<typename Member, int Ndim>
vector<hpx::id_type,tree<Member,Ndim>::Nbranch> tree<Member,Ndim>::child_gids;

template<typename Member, int Ndim>
hpx::id_type tree<Member,Ndim>::this_gid;

template<typename Member, int Ndim>
tree<Member,Ndim>* tree<Member,Ndim>::this_ptr;

}




#endif /* TREE_HPP_ */
