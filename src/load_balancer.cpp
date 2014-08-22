/*
 * load_balancer.cpp
 *
 *  Created on: May 31, 2014
 *      Author: dmarce1
 */

#include <hpx/config.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/list.hpp>
#include "load_balancer.hpp"

namespace xtree {

load_balancer::load_balancer() {
	assert(false);
}

load_balancer::load_balancer(component_type* back_ptr) :
		base_type(back_ptr), mins_begin(INT_MAX, hpx::invalid_id), semaphore(1) {
	static bool initialized = false;
	assert(!initialized);
//	printf( "A\n");
	initialized = true;
	int mask, neighbor_id, my_id, id_cnt;
	auto localities = hpx::find_all_localities();
	std::vector<hpx::future<hpx::id_type>> futures;
	my_load = 0;
	my_id = hpx::get_locality_id();
//	printf( "B\n");
	id_cnt = localities.size();
	hpx::register_id_with_basename(name, base_type::get_gid(), my_id).get();
	std::list<int> tmp;
	tmp.push_front(my_id);
	for (mask = 1; mask > 0; mask <<= 1) {
		neighbor_id = mask ^ my_id;
		if (neighbor_id < id_cnt) {
			tmp.push_front(neighbor_id);
		}
	}
//	printf( "C\n");
	std::vector<int> tmp2(tmp.size());
	std::copy(tmp.begin(), tmp.end(), tmp2.begin());
	std::sort(tmp2.begin(), tmp2.end());
	neighbors.resize(tmp2.size());
	futures.resize(neighbors.size());
	for (size_t i = 0; i < neighbors.size(); i++) {
	//	printf( "%i\n", tmp2[i]);
		futures[i] = hpx::find_id_from_basename(name, tmp2[i]);
	}
	//printf( "D\n");
	for (size_t i = 0; i < neighbors.size(); i++) {
//		printf( "----D\n");
		neighbors[i] = futures[i].get();
	}
//	printf( "E\n");

}

load_balancer::~load_balancer() {
	// TODO Auto-generated destructor stub
}

void load_balancer::decrement_load() {
	semaphore.wait();
	my_load--;
	semaphore.signal();
}

int load_balancer::unlock_servlet(bool inc_cnt) {

	if (inc_cnt) {
		my_load++;
	}
	semaphore.signal();
	return hpx::get_locality_id();
}

hpx::future<std::pair<int, hpx::id_type>> load_balancer::lock_servlet(std::pair<int, hpx::id_type> best_mins, std::list<hpx::id_type> remaining) {
	hpx::future<std::pair<int, hpx::id_type>> fut;

	semaphore.wait();
	if (my_load < best_mins.first) {
		if (best_mins.second != hpx::invalid_id) {
			hpx::apply<action_unlock_servlet>(best_mins.second, false);
		}
		best_mins.first = my_load;
		best_mins.second = base_type::get_gid();
	} else {
		semaphore.signal();
	}

	if (remaining.size() > 0) {
		const auto id = remaining.front();
		remaining.pop_front();
		fut = hpx::async<action_lock_servlet>(id, best_mins, std::move(remaining));
	} else {
		fut = hpx::make_ready_future(best_mins);
	}

	return fut;
}

hpx::future<int> load_balancer::increment_load() {

	std::pair<int, hpx::id_type> mins;
	mins.first = INT_MAX;
	mins.second = hpx::invalid_id;
	std::list<hpx::id_type> remaining;
	for (size_t i = 1; i < neighbors.size(); i++) {
		remaining.push_back(neighbors[i]);
	}
	auto fut1 = hpx::async<action_lock_servlet>(neighbors[0], mins_begin, remaining);
	auto fut2 = fut1.then(hpx::util::unwrapped([](std::pair<int, hpx::id_type> mins) {
		assert(mins.second != hpx::invalid_id);
		return hpx::async<action_unlock_servlet>(mins.second, true);
	}));
	return fut2;
}

load_balancer* load_balancer::get_ptr() {
	return this;
}

int load_balancer::get_load() {
	int i;
	semaphore.wait();
	i = my_load;
	semaphore.signal();
	return i;
}
} /* namespace xtree */
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::load_balancer>, load_balancer);

HPX_REGISTER_ACTION(lock_servlet_action);
HPX_REGISTER_ACTION(unlock_servlet_action);
HPX_REGISTER_ACTION(get_ptr_action);


