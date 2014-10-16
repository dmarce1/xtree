/*
 * load_balancer.cpp
 *
 *  Created on: May 31, 2014
 *      Author: dmarce1
 */

#include <hpx/config.hpp>
#include <boost/serialization/list.hpp>
#include "load_balancer.hpp"

namespace xtree {

load_balancer::load_balancer() {
	assert(false);
}

load_balancer::load_balancer(component_type* back_ptr) :
		base_type(back_ptr) {
	const char name[] = "load_balancer";
	hpx::register_id_with_basename(name, base_type::get_gid(), hpx::get_locality_id()).get();
	auto hfut = hpx::find_id_from_basename(name, 0);

	if (hpx::get_locality_id() == 0) {
		auto localities = hpx::find_all_localities();
		procs.resize(localities.size());
	}

	home = hfut.get();

}

load_balancer::~load_balancer() {
	printf("Delete\n");
	hpx::unregister_id_with_basename(name, hpx::get_locality_id()).get();
}

void load_balancer::decrement_load(const location<Ndim>& loc) {
	hpx::async < decrement_server_action > (home, hpx::get_locality_id(), loc).get();
}

hpx::future<int> load_balancer::increment_load(const location<Ndim>& l) {
	return hpx::async < increment_server_action > (home, l);
}

load_balancer* load_balancer::get_ptr() {
	return this;
}

void load_balancer::print() {

	printf("\n----------Load Balance---------------\n");
	for (int i = 0; i < procs.size(); i++) {
		printf("%6i %6li\n", i, procs[i].list.size());
	}
	printf("\n");
}

int load_balancer::increment_server(const location<Ndim>& this_loc) {
	boost::lock_guard < hpx::lcos::local::spinlock > this_lock(lock);
	double load_score, proximity_score, neighbor_score;
	std::list<int> eligible_procs;
	int lowest_load = INT_MAX;
	int highest_load = INT_MIN;
	double best_score = -DBL_MAX;
	double score;
	int next_proc;
	for (int i = 0; i < procs.size(); i++) {
		lowest_load = std::min(lowest_load, int(procs[i].list.size()));
		highest_load = std::max(highest_load, int(procs[i].list.size()));
	}

	int window = 1 + lowest_load / 2;

	for (int i = 0; i < procs.size(); i++) {
		if (procs[i].list.size() - lowest_load <= window) {
			eligible_procs.push_back(i);
		}
	}
	location<Ndim> search_loc;
	for (auto i = eligible_procs.begin(); i != eligible_procs.end(); ++i) {
		int neighbor_cnt = 0;
		for (int j = 0; j < Ndim; j++) {
			dir_type<Ndim> dir;
			dir.set_zero();
			dir[j] = +1;
			if (!this_loc.is_phys_bnd(dir)) {
				search_loc = this_loc.get_neighbor(dir);
				if (procs[*i].list.find(std::move(search_loc)) != procs[*i].list.end()) {
					neighbor_cnt++;
				}
			}
			dir[j] = -1;
			if (!this_loc.is_phys_bnd(dir)) {
				search_loc = this_loc.get_neighbor(dir);
				if (procs[*i].list.find(std::move(search_loc)) != procs[*i].list.end()) {
					neighbor_cnt++;
				}
			}
		}
		if (this_loc.get_level() > 0) {
			search_loc = this_loc.get_parent();
			if (procs[*i].list.find(std::move(search_loc)) != procs[*i].list.end()) {
				neighbor_cnt++;
			}
		}
		neighbor_score = double(neighbor_cnt) / double(2 * Ndim + 1);

		double d = 0.0;
		auto tmp = this_loc.get_position();
		for (int j = 0; j < Ndim; j++) {
			d += pow(tmp[j] - procs[*i].avg_loc[j], 2);
		}
		d = sqrt(d / double(Ndim));
		proximity_score = 1.0 - d;
//		printf( "%e\n", proximity_score);

		load_score = double(highest_load - procs[*i].list.size()+1) / double(highest_load - lowest_load+1);
		score = load_score + sqrt(10.0)*proximity_score + 10.0*neighbor_score;
		if (score > best_score) {
			best_score = score;
			next_proc = *i;
		}
	}
	if (rand() % 100 == 0) {
		print();
	}
	auto& e = procs[next_proc];
	if (e.list.size() == 0) {
		for (int i = 0; i < Ndim; i++) {
			e.avg_loc[i] = 0.0;
		}
	} else {
		for (int i = 0; i < Ndim; i++) {
			e.avg_loc[i] *= e.list.size();
		}
	}
	for (int i = 0; i < Ndim; i++) {
		e.avg_loc[i] += this_loc.get_position()[i];
	}
	procs[next_proc].list.insert(this_loc);
	for (int i = 0; i < Ndim; i++) {
		e.avg_loc[i] /= e.list.size();
	}
	return next_proc;
}

void load_balancer::decrement_server(int proc, const location<Ndim>& loc) {
	boost::lock_guard < hpx::lcos::local::spinlock > this_lock(lock);
	auto i = procs[proc].list.find(loc);
	assert(i != procs[proc].list.end());
	auto& e = procs[proc];
	for (int i = 0; i < Ndim; i++) {
		e.avg_loc[i] *= e.list.size();
	}
	for (int i = 0; i < Ndim; i++) {
		e.avg_loc[i] -= loc.get_position()[i];
	}
	procs[proc].list.erase(std::move(i));
	for (int i = 0; i < Ndim; i++) {
		e.avg_loc[i] /= e.list.size();
	}
}

int load_balancer::get_load() {
	int i;
	return i;
}
} /* namespace xtree */
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::load_balancer>, load_balancer);

HPX_REGISTER_ACTION (increment_server_action);
HPX_REGISTER_ACTION (decrement_server_action);
//HPX_REGISTER_ACTION(unlock_servlet_action);
HPX_REGISTER_ACTION (get_ptr_action);

