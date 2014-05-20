#include "xtree.hpp"
#include "vector.hpp"

namespace xtree {
namespace server {

static std::vector<hpx::id_type> neighbors;
static hpx::lcos::local::counting_semaphore semaphore(1);
static hpx::lcos::local::counting_semaphore lock(1);
static int my_load = 0;

static std::pair<int, hpx::id_type> lock_servlet(std::pair<int, hpx::id_type>);
static void unlock_servlet(bool);
static void locality_server_initialize();

using action_lock_servlet = HPX_MAKE_ACTION(xtree::server::lock_servlet)::type;
using action_unlock_servlet = HPX_MAKE_ACTION(xtree::server::unlock_servlet)::type;

int get_load() {
	int i;
	semaphore.wait();
	i = my_load;
	semaphore.signal();
	return i;
}

hpx::future<hpx::id_type> increment_load() {
	locality_server_initialize();
	std::pair<int, hpx::id_type> mins;
	mins.first = INT_MAX;
	mins.second = hpx::invalid_id;
	auto fut = hpx::make_ready_future(mins);
	for (size_t i = 0; i < neighbors.size(); i++) {
		fut = fut.then(hpx::util::unwrapped([i](std::pair<int,hpx::id_type> mins) {
			return hpx::async<action_lock_servlet>(neighbors[i], mins);
		}));
	}
	return fut.then(hpx::util::unwrapped([](std::pair<int,hpx::id_type> mins) {
		hpx::apply<action_unlock_servlet>(mins.second, true);
		return mins.second;
	}));
}

void decrement_load() {
	semaphore.wait();
	my_load--;
	semaphore.signal();
}

static void locality_server_initialize() {
	static bool initialized = false;
	if (!initialized) {
		lock.wait();
		if (!initialized) {
			std::list<int> tmp;
			int mask, neighbor_id;
			const int my_id = hpx::get_locality_id();
			std::vector<hpx::id_type> localities = hpx::find_all_localities();
			const int id_cnt = localities.size();
			tmp.push_back(my_id);
			for (mask = 1; mask > 0; mask <<= 1) {
				neighbor_id = mask ^ my_id;
				if (neighbor_id < id_cnt) {
					tmp.push_back(neighbor_id);
				}
			}
			std::vector<int> tmp2 = std::vector<int>(tmp.begin(), tmp.end());
			std::sort(tmp2.begin(), tmp2.end());
			neighbors.resize(tmp2.size());
			for (size_t i = 0; i < neighbors.size(); i++) {
				neighbors[i] = localities[tmp2[i]];
			}
			initialized = true;
		}
		lock.signal();
	}
}

static std::pair<int, hpx::id_type> lock_servlet(std::pair<int, hpx::id_type> mins) {
	locality_server_initialize();
	semaphore.wait();
	if (my_load < mins.first) {
		mins.first = my_load;
		if (mins.second != hpx::invalid_id) {
			hpx::apply<action_unlock_servlet>(mins.second, false);
		}
		mins.second = hpx::find_here();
	} else {
		semaphore.signal();
	}
	return mins;
}

static void unlock_servlet(bool inc_cnt) {
	locality_server_initialize();
	if (inc_cnt) {
		my_load++;
	}
	semaphore.signal();
}

}
}

typedef xtree::server::action_lock_servlet xtree_server_action_lock_servlet;
typedef xtree::server::action_unlock_servlet xtree_server_action_unlock_servlet;
HPX_REGISTER_PLAIN_ACTION(xtree_server_action_lock_servlet);
HPX_REGISTER_PLAIN_ACTION(xtree_server_action_unlock_servlet);
