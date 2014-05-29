/*
 * xtree.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef XTREE_HPP_
#define XTREE_HPP_

#include "fwd.hpp"

#define XTREE_MAKE_ACTION( a, b ) 																					\
using a = typename hpx::actions::make_action<decltype(&b),&b>::type;	 											\
/**/

#define XTREE_INSTANTIATE( DERIVED_CLASS, ... )																		\
namespace xtree {																									\
	typedef node<DERIVED_CLASS, __VA_ARGS__> base_node_type;   														\
	typedef tree<DERIVED_CLASS, __VA_ARGS__> tree_type;   															\
}																													\
typedef xtree::tree<DERIVED_CLASS,__VA_ARGS__>::action_get_new_node action_get_new_node_global;	\
typedef xtree::tree<DERIVED_CLASS,__VA_ARGS__>::action_get_max_level action_get_max_level_global;	\
typedef xtree::tree<DERIVED_CLASS,__VA_ARGS__>::action_get_terminal_future action_get_terminal_future_global;	\
typedef xtree::tree<DERIVED_CLASS,__VA_ARGS__>::action_lock_servlet action_lock_servlet_global;	\
typedef xtree::tree<DERIVED_CLASS,__VA_ARGS__>::action_unlock_servlet action_unlock_servlet_global;	\
HPX_REGISTER_PLAIN_ACTION(action_get_new_node_global);	\
HPX_REGISTER_PLAIN_ACTION(action_get_max_level_global);	\
HPX_REGISTER_PLAIN_ACTION(action_get_terminal_future_global);	\
HPX_REGISTER_PLAIN_ACTION(action_lock_servlet_global);	\
HPX_REGISTER_PLAIN_ACTION(action_unlock_servlet_global);	\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::tree_type>, tree_type);	\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::base_node_type>, base_node_type);	\
/**/


#include <hpx/hpx_init.hpp>
#include <hpx/lcos/local/counting_semaphore.hpp>
#include <hpx/include/components.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/local/mutex.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/util/unwrapped.hpp>

#include <boost/mpl/int.hpp>
#include <boost/serialization/list.hpp>

#include <silo.h>

#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <unordered_set>
#include <utility>
#include <vector>


#include "util.hpp"
#include "vector.hpp"
#include "grid_base.hpp"
#include "grid.hpp"
#include "bgrid.hpp"
#include "indexer.hpp"
#include "location.hpp"
#include "node.hpp"
#include "tree.hpp"


#endif /* XTREE_HPP_ */
