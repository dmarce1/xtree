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

	if( hpx::get_locality_id() == 0 ) {
		auto localities = hpx::find_all_localities();
		procs.resize(localities.size());
	}


	home = hfut.get();


}

load_balancer::~load_balancer() {
	printf( "Delete\n");
	hpx::unregister_id_with_basename(name, hpx::get_locality_id()).get();
}

void load_balancer::decrement_load(const location<Ndim>& loc) {
	hpx::async<decrement_server_action>(home,hpx::get_locality_id(), loc).get();
}

hpx::future<int> load_balancer::increment_load(const location<Ndim>& l) {
	return hpx::async<increment_server_action>(home, l);
}

load_balancer* load_balancer::get_ptr() {
	return this;
}

void load_balancer::print() {

		printf( "\n----------Load Balance---------------\n");
		for( int i = 0; i < procs.size(); i++) {
			printf( "%6i %6li\n",i, procs[i].size());
		}
		printf( "\n");
}



int load_balancer::increment_server(const location<Ndim>& this_loc) {
	boost::lock_guard<hpx::lcos::local::spinlock> this_lock(lock);

	const int window = 8;
	std::list<int> eligible_procs;
	int lowest_load = INT_MAX;
	int most_neighbors = -1;
	int next_proc;
	for( int i =0; i < procs.size(); i++) {
		lowest_load = std::min(lowest_load, int(procs[i].size()));
	}
	for( int i =0; i < procs.size(); i++) {
		if( procs[i].size() - lowest_load <= window ) {
			eligible_procs.push_back(i);
		}
	}
	for( auto i = eligible_procs.begin(); i != eligible_procs.end(); ++i) {
		int neighbor_cnt = 0;
		for( dir_type<Ndim> dir; !dir.is_end(); dir++) {
			location<Ndim> search_loc = this_loc.get_neighbor(dir);
			if( procs[*i].find(std::move(search_loc)) != procs[*i].end()) {
				neighbor_cnt++;
			}
		}
		for( child_index_type<Ndim> ci; !ci.is_end(); ci++) {
			location<Ndim> search_loc = this_loc.get_child(ci);
			if( procs[*i].find(std::move(search_loc)) != procs[*i].end()) {
				neighbor_cnt++;
			}
		}
		location<Ndim> search_loc = this_loc.get_parent();
		if( procs[*i].find(std::move(search_loc)) != procs[*i].end()) {
			neighbor_cnt++;
		}
		if( neighbor_cnt > most_neighbors ) {
			most_neighbors = neighbor_cnt;
			next_proc = *i;
		}
	}
	if( rand() % 100 == 0 ) {
		print();
	}
	procs[next_proc].insert(this_loc);
	return next_proc;
}

void load_balancer::decrement_server(int proc,const location<Ndim>& loc) {
	boost::lock_guard<hpx::lcos::local::spinlock> this_lock(lock);
	auto i = procs[proc].find(loc);
	assert( i != procs[proc].end());
	procs[proc].erase(std::move(i));
}


int load_balancer::get_load() {
	int i;
	return i;
}
} /* namespace xtree */
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::load_balancer>, load_balancer);

HPX_REGISTER_ACTION(increment_server_action);
HPX_REGISTER_ACTION(decrement_server_action);
//HPX_REGISTER_ACTION(unlock_servlet_action);
HPX_REGISTER_ACTION(get_ptr_action);


