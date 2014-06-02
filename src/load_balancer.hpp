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

namespace xtree {

class load_balancer: public hpx::components::managed_component_base<load_balancer, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr> {
public:
	const char* name = "load_balancer";
	using base_type = hpx::components::managed_component_base<load_balancer, hpx::components::detail::this_type, hpx::traits::construct_with_back_ptr>;
	using component_type = hpx::components::managed_component<load_balancer>;
private:
	const std::pair<int, hpx::id_type> mins_begin;
	hpx::lcos::local::counting_semaphore semaphore;
	int my_load;
	std::vector<hpx::id_type> neighbors;
public:
	load_balancer();
	load_balancer(component_type*);
	virtual ~load_balancer();
	template<typename Archive>
	void serialize(Archive& ar, const int v) {
		assert(false);
	}
	hpx::future<hpx::id_type> increment_load();
	void decrement_load();
	hpx::id_type unlock_servlet(bool inc_cnt);
	hpx::future<std::pair<int, hpx::id_type>> lock_servlet(std::pair<int, hpx::id_type> best_mins, std::list<hpx::id_type> remaining);
	int get_load();
	load_balancer* get_ptr();
	using action_lock_servlet = typename hpx::actions::make_action<decltype(&load_balancer::lock_servlet),&load_balancer::lock_servlet>::type;
	using action_unlock_servlet = typename hpx::actions::make_action<decltype(&load_balancer::unlock_servlet),&load_balancer::unlock_servlet>::type;
	using action_get_ptr = typename hpx::actions::make_action<decltype(&load_balancer::get_ptr),&load_balancer::get_ptr>::type;

};

} /* namespace xtree */

#endif /* LOAD_BALANCER_H_ */
