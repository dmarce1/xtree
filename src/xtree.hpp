/*
 * xtree.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef XTREE_HPP_
#define XTREE_HPP_

#include "int_seq.hpp"
#include "load_balancer.hpp"
#include "node.hpp"

namespace xtree {
template<typename Derived, typename Dims>
class tree: public hpx::components::managed_component_base<tree<Derived, Dims>> {
public:
	static constexpr int Nbranch = 2;
	static constexpr int Ndim = Dims::dim();
	static constexpr int Nneighbor = pow_<3, Ndim>::get();
private:
	using node_type = node<Derived, Dims>;
	std::set<Derived*> nodes;
	vector<hpx::id_type, Nbranch> child_ids;
public:
	template<typename Arc>
	void serialize(Arc& arc, const int v) {
		assert(false);
	}
	tree() {
		std::vector<hpx::future<hpx::id_type>> futures(Nbranch);
		auto localities = hpx::find_all_localities();
		int myloc = hpx::get_locality_id();
		if (myloc == 0) {
			load_balancer::initialize().get();
		}
		for (int i = 0; i < Nbranch; i++) {
			int j = myloc * Nbranch + i + 1;
			if (j < localities.size()) {
				futures[j] = hpx::new_<tree<Derived, Dims>>(localities[j]);
			} else {
				futures[j] = hpx::make_ready_future(hpx::invalid_id);
			}
		}
		for (int i = 0; i < Nbranch; i++) {
			child_ids[i] = futures[i].get();
		}
		if (myloc == 0) {
			load_balancer::increment_load();
			new_node(location<Ndim>(), hpx::invalid_id, vector<hpx::id_type, Nneighbor>(hpx::invalid_id));
		}
	}
	virtual ~tree() {
	}
	hpx::id_type new_node(const location<Ndim>& _loc, hpx::id_type _parent_id, vector<hpx::id_type, Nneighbor> _neighbors) {
		hpx::id_type newnode = hpx::new_<node_type>(hpx::find_here(), this, _loc, _parent_id, _neighbors).get();
		nodes.insert(hpx::async<typename node<Derived, Dims>::action_get_this>(newnode).get());
		return newnode;
	}
	void delete_node(Derived* ptr) {
		nodes.erase(nodes.find(ptr));
	}

	template<enum node_type::op_type Op, typename T, node<Derived, Dims>::get_type<T> Get, node<Derived, Dims>::set_type<T> Set>
	void begin_execute() {
		int lev;
		switch (Op) {
		case node<Derived, Dims>::ASCEND:
			lev = 0;
			while (execute<Op, T, Get, Set>(lev).get()) {
				lev++;
			}
			break;
		case node<Derived, Dims>::DESCEND:
			lev = get_max_level().get();
			while (lev >= 0) {
				execute<Op, T, Get, Set>(lev).get();
				lev--;
			}
			break;
		case node<Derived, Dims>::EXCHANGE:
			execute<Op, T, Get, Set>(-1).get();
			break;
		}

	}

	template<enum node_type::op_type Op, typename T, node<Derived, Dims>::get_type<T> Get, node<Derived, Dims>::set_type<T> Set>
	hpx::future<bool> execute(int level) {
		bool rc = false;
		std::vector<hpx::future<bool>> futures(Nbranch + nodes.size());
		for (int i = 0; i < Nbranch; i++) {
			if (child_ids[i] != hpx::invalid_id) {
				futures[i] = hpx::async<typename operations<Op, T, Get, Set>::action_execute>(child_ids[i], level);
			} else {
				futures[i] = hpx::make_ready_future(true);
			}
		}
		for (auto i = nodes.begin(); i != nodes.end(); i++) {
			if (level == -1 || ((*i)->get_level() == level)) {
				rc = true;
				(*i)->template setup_op_dataflow<Op, T, Get, Set>();
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

	hpx::future<int> get_max_level() {
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
	using action_get_max_level = typename hpx::actions::make_action<decltype(&tree::get_max_level ), &tree::get_max_level >::type;

	template<enum node_type::op_type Op, typename T, node<Derived, Dims>::get_type<T> Get, node<Derived, Dims>::set_type<T> Set>
	struct operations {
		using action_execute = typename hpx::actions::make_action<decltype(&tree::execute<Op,T,Get,Set>), &tree::execute<Op,T,Get,Set>>::type;

	};
};
}

#define XTREE_INSTANTIATE( DERIVED_CLASS, ... )																		\
namespace xtree {																									\
	typedef node<DERIVED_CLASS, int_seq<__VA_ARGS__>> base_node_type;   											\
	typedef tree<DERIVED_CLASS, int_seq<__VA_ARGS__>> tree_type;   											\
}																													\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::base_node_type>, base_node_type);	\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::tree_type>, tree_type);	\
/**/

#endif /* XTREE_HPP_ */
