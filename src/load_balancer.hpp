/*
 * load_server.hpp
 *
 *  Created on: May 20, 2014
 *      Author: dmarce1
 */

#ifndef load_server_HPP_
#define load_server_HPP_

#include <hpx/include/components.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/lcos/local/counting_semaphore.hpp>


namespace xtree {
class load_balancer {
private:
	static const std::pair<int, hpx::id_type> mins_begin;
	static std::vector<hpx::id_type> neighbors;
	static hpx::lcos::local::counting_semaphore semaphore;
	static int my_load;

	static hpx::future<std::pair<int, hpx::id_type>> lock_servlet(std::pair<int, hpx::id_type>, std::list<hpx::id_type> remaining);
	static hpx::id_type unlock_servlet(bool);
	static hpx::future<void> initialize_(std::list<hpx::id_type>);

public:
	static hpx::future<void> initialize();
	static hpx::future<hpx::id_type> increment_load();
	static int get_load();
	static void decrement_load();

	using action_lock_servlet = HPX_MAKE_ACTION(xtree::load_balancer::lock_servlet)::type;

	using action_unlock_servlet = HPX_MAKE_ACTION(xtree::load_balancer::unlock_servlet)::type;

	using action_initialize_ = HPX_MAKE_ACTION(xtree::load_balancer::initialize_)::type;

};
}

#endif /* load_server_HPP_ */
