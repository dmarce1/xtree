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
		for( std::size_t i = 0; i < localities.size(); i++) {
			procs[i].load = 0;
			procs[i].loc = std::vector<double>(Ndim,0.0);
		}
	}


	home = hfut.get();


}

load_balancer::~load_balancer() {
	hpx::unregister_id_with_basename(name, hpx::get_locality_id()).get();
}

void load_balancer::decrement_load(const std::array<double, Ndim>& loc) {
	hpx::apply<decrement_server_action>(home,hpx::get_locality_id(), loc);
}

hpx::future<int> load_balancer::increment_load(const std::array<double, Ndim>& l) {
	return hpx::async<increment_server_action>(home, l);
}

load_balancer* load_balancer::get_ptr() {
	return this;
}

void load_balancer::print() {

		printf( "\n----------Load Balance---------------\n");
		for( int i = 0; i < procs.size(); i++) {
			printf( "%6i %6i %e %e %e\n",i, procs[i].load, procs[i].loc[0], procs[i].loc[1], procs[i].loc[2]);
		}
		printf( "\n");
}



int load_balancer::increment_server(const std::array<double, Ndim>& pos0) {
	boost::lock_guard<hpx::lcos::local::spinlock> this_lock(lock);
	int next_proc;

	double best_score = 1.0e+99;
	int best_load = INT_MAX;

	for( int i = 0; i < procs.size(); i++ ) {
		best_load = std::min(best_load,procs[i].load);
	}

	for( int i = 0; i < procs.size(); i++ ) {
		double score = 0.0;
		for( int j = 0; j < Ndim; j++) {
			score += pow(pos0[j] - procs[i].loc[j],2);
		}
		score = sqrt(score) / sqrt(Ndim);
		score += 4.0*((procs[i].load + 1.0) / (best_load + 1.0) - 1.0);
		if( score < best_score ) {
			next_proc = i;
			best_score = score;
		}
	}

	for( int i = 0; i < Ndim; i++) {
		procs[next_proc].loc[i] *= procs[next_proc].load;
		procs[next_proc].loc[i] += pos0[i];
	}
	procs[next_proc].load++;
	for( int i = 0; i < Ndim; i++) {
		procs[next_proc].loc[i] /= procs[next_proc].load;
	}

	if( rand() % 100 == 0 ) {
		print();
	}
	return next_proc;
}

void load_balancer::decrement_server(int proc,const std::array<double, Ndim>& pos0) {
	boost::lock_guard<hpx::lcos::local::spinlock> this_lock(lock);

	for( int i = 0; i < Ndim; i++) {
		procs[proc].loc[i] *= procs[proc].load;
		procs[proc].loc[i] -= pos0[i];
	}
	procs[proc].load--;
	if( procs[proc].load != 0 ) {
		for( int i = 0; i < Ndim; i++) {
			procs[proc].loc[i] /= procs[proc].load;
		}
	}
	for( int i = 0; i < Ndim; i++) {
		procs[proc].loc[i] = 0.0;
	}


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


