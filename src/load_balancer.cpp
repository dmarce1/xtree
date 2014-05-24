#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/util/unwrapped.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include "load_balancer.hpp"

namespace xtree {

const std::pair<int, hpx::id_type> load_balancer::mins_begin = std::pair<int, hpx::id_type>(INT_MAX, hpx::invalid_id);
std::vector<hpx::id_type> load_balancer::neighbors;
hpx::lcos::local::counting_semaphore load_balancer::semaphore(1);
int load_balancer::my_load = 0;

hpx::future<void> load_balancer::initialize() {
	auto localities = hpx::find_all_localities();
	std::list<hpx::id_type> list(localities.begin(), localities.end());
	return hpx::async<action_initialize_>(list.front(), list);
}

hpx::future<std::pair<int, hpx::id_type>> load_balancer::lock_servlet(std::pair<int, hpx::id_type> best_mins, std::list<hpx::id_type> remaining) {
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

hpx::id_type load_balancer::unlock_servlet(bool inc_cnt) {

	if (inc_cnt) {
		my_load++;
	}
	semaphore.signal();
	return hpx::find_here();
}

int load_balancer::get_load() {
	int i;
	semaphore.wait();
	i = my_load;
	semaphore.signal();
	return i;
}

hpx::future<hpx::id_type> load_balancer::increment_load() {

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

void load_balancer::decrement_load() {
	semaphore.wait();
	my_load--;
	semaphore.signal();
}

hpx::future<void> load_balancer::initialize_(std::list<hpx::id_type> remaining) {
	std::list<hpx::id_type> lists[2];
	std::list<int> tmp;
	std::vector<hpx::future<void>> futs(2);
	int mask, neighbor_id, my_id, id_cnt;

	remaining.pop_front();
	for (int j = 0; !remaining.empty(); j ^= 1) {
		lists[j].push_back(remaining.front());
		remaining.pop_front();
	}
	for (int i = 0; i < 2; i++) {
		if (lists[i].size() > 0) {
			futs[i] = hpx::async<action_initialize_>(lists[i].front(), lists[i]);
		} else {
			futs[i] = hpx::make_ready_future();
		}
	}
	my_id = hpx::get_locality_id();
	const auto localities = hpx::find_all_localities();
	id_cnt = localities.size();
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
		neighbors[i] = localities[tmp2[i]];
	}
	return when_all(futs);
}

}

typedef xtree::load_balancer::action_lock_servlet xtree_load_balancer_action_lock_servlet;
typedef xtree::load_balancer::action_unlock_servlet xtree_load_balancer_action_unlock_servlet;
typedef xtree::load_balancer::action_initialize_ xtree_load_balancer_action_initialize_;
HPX_REGISTER_PLAIN_ACTION(xtree_load_balancer_action_lock_servlet);
HPX_REGISTER_PLAIN_ACTION(xtree_load_balancer_action_unlock_servlet);
HPX_REGISTER_PLAIN_ACTION(xtree_load_balancer_action_initialize_);
