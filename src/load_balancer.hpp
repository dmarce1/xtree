/*
 * load_balancer.h
 *
 *  Created on: May 31, 2014
 *      Author: dmarce1
 */

#ifndef LOAD_BALANCER_H_
#define LOAD_BALANCER_H_

#include <hpx/include/components.hpp>
#include <hpx/lcos/local/counting_semaphore.hpp>
#include <hpx/lcos/local/mutex.hpp>
#include "location.hpp"
#include <unordered_set>

namespace xtree {

class load_balancer: public hpx::components::managed_component_base<load_balancer, hpx::components::detail::this_type,
		hpx::traits::construct_with_back_ptr> {
	static constexpr int Ndim = 3;
public:
	const char* name = "load_balancer";
	using base_type = hpx::components::managed_component_base<load_balancer, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr>;
	using component_type = hpx::components::managed_component<load_balancer>;
private:
	hpx::id_type home;
	mutable hpx::lcos::local::spinlock lock;
	struct  entry_type{
		std::unordered_set<location<Ndim>> list;
		std::array<double,Ndim> avg_loc;
	} ;
	std::vector<entry_type> procs;
public:
	void print();
	load_balancer();
	load_balancer(component_type*);
	virtual ~load_balancer();
	template<typename Archive>
	void serialize(Archive& ar, const int v) {
		assert(false);
	}

	hpx::future<int> increment_load(const location<Ndim>&);
	void decrement_load(const location<Ndim>&);
	int increment_server(const location<Ndim>&);
	void decrement_server(int, const location<Ndim> &);
	int get_load();
	load_balancer* get_ptr();
	using action_increment_server = typename hpx::actions::make_action<decltype(&load_balancer::increment_server),&load_balancer::increment_server>::type;
	using action_decrement_server = typename hpx::actions::make_action<decltype(&load_balancer::decrement_server),&load_balancer::decrement_server>::type;
	using action_get_ptr = typename hpx::actions::make_action<decltype(&load_balancer::get_ptr),&load_balancer::get_ptr>::type;

};

} /* namespace xtree */

typedef xtree::load_balancer::action_increment_server increment_server_action;
typedef xtree::load_balancer::action_decrement_server decrement_server_action;
typedef xtree::load_balancer::action_get_ptr get_ptr_action;

HPX_REGISTER_ACTION_DECLARATION (increment_server_action);
HPX_REGISTER_ACTION_DECLARATION (decrement_server_action);
//HPX_REGISTER_ACTION_DECLARATION(unlock_servlet_action);
HPX_REGISTER_ACTION_DECLARATION (get_ptr_action);

#endif /* LOAD_BALANCER_H_ */
