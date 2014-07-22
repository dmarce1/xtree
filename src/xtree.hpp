/*
 * xtree.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef XTREE_HPP_
#define XTREE_HPP_

#include <hpx/hpx_init.hpp>
#include "fwd.hpp"
#include "container_math.hpp"
#include <hpx/runtime/components/derived_component_factory.hpp>

#define XTREE_INSTANTIATE( DERIVED_CLASS, ... )																		\
namespace xtree {																									\
	typedef silo_output<__VA_ARGS__> silo_output_type;   														\
	typedef node<DERIVED_CLASS, __VA_ARGS__> node_type;   														\
	typedef tree<DERIVED_CLASS, __VA_ARGS__> tree_type;   															\
	typedef  hpx::components::managed_component<xtree::node_type> node_component_type; \
	typedef DERIVED_CLASS derived_type;\
}																													\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::silo_output_type>, silo_output_type);	\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::tree_type>, tree_type);	\
HPX_DEFINE_GET_COMPONENT_TYPE(xtree::node_type);\
HPX_REGISTER_DERIVED_COMPONENT_FACTORY(hpx::components::managed_component<DERIVED_CLASS>, DERIVED_CLASS, "node_type");	\
typedef DERIVED_CLASS::action_operations_end action0_operations_end; \
HPX_REGISTER_ACTION_DECLARATION(action0_operations_end); \
HPX_REGISTER_ACTION(action0_operations_end); \
typedef DERIVED_CLASS::action_initialize action0_initialize; \
HPX_REGISTER_ACTION_DECLARATION(action0_initialize); \
HPX_REGISTER_ACTION(action0_initialize); \
typedef DERIVED_CLASS::action_debranch action0_debranch; \
HPX_REGISTER_ACTION_DECLARATION(action0_debranch); \
HPX_REGISTER_ACTION(action0_debranch); \
typedef DERIVED_CLASS::action_get_this action0_get_this; \
HPX_REGISTER_ACTION_DECLARATION(action0_get_this); \
HPX_REGISTER_ACTION(action0_get_this); \
typedef DERIVED_CLASS::action_find_neighbors action0_find_neighbors; \
HPX_REGISTER_ACTION_DECLARATION(action0_find_neighbors); \
HPX_REGISTER_ACTION(action0_find_neighbors); \
typedef DERIVED_CLASS::action_notify_branch action0_notify_branch; \
HPX_REGISTER_ACTION_DECLARATION(action0_notify_branch); \
HPX_REGISTER_ACTION(action0_notify_branch); \
typedef DERIVED_CLASS::action_notify_debranch action0_notify_debranch; \
HPX_REGISTER_ACTION_DECLARATION(action0_notify_debranch); \
HPX_REGISTER_ACTION(action0_notify_debranch); \
typedef DERIVED_CLASS::action_note_sibs action0_note_sibs0; \
HPX_REGISTER_ACTION_DECLARATION(action0_note_sibs0); \
HPX_REGISTER_ACTION(action0_note_sibs0); \
typedef DERIVED_CLASS::action_notify_of_neighbor action0_notify_of_neighbor; \
HPX_REGISTER_ACTION_DECLARATION(action0_notify_of_neighbor); \
HPX_REGISTER_ACTION(action0_notify_of_neighbor); \
typedef xtree::tree_type::action_get_new_node action_get_new_node; \
HPX_REGISTER_ACTION (action_get_new_node); \
typedef xtree::tree_type::action_get_this action_get_this;  \
HPX_REGISTER_ACTION (action_get_this); \
typedef xtree::tree_type::action_place_root action_place_root;  \
HPX_REGISTER_ACTION (action_place_root); \
typedef xtree::tree_type::action_output action_output; \
HPX_REGISTER_ACTION (action_output); \
typedef xtree::silo_output_type::action_send_zones_to_silo action_send_zones_to_silo; \
HPX_REGISTER_ACTION( action_send_zones_to_silo );


/**/

#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/util/unwrapped.hpp>

#include <boost/serialization/valarray.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/mpl/int.hpp>

#include <silo.h>

#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <unordered_set>
#include <utility>
#include <vector>

namespace boost {
namespace serialization {

template<class Archive, class T, size_t N>
void serialize(Archive & ar, std::array<T, N> & a, const unsigned int version) {
	ar & boost::serialization::make_array(a.data(), a.size());
}

/*template<class Archive, class T>
void serialize(Archive & ar, std::vector<T> & a, const unsigned int version) {
	ar & boost::serialization::make_array(a.data(), a.size());
}

template<class Archive, class T>
void serialize(Archive & ar, const std::vector<T> & a, const unsigned int version) {
	ar << boost::serialization::make_array(a.data(), a.size());
}
*/
} // namespace serialization
} // namespace boost
#include "load_balancer.hpp"
#include "util.hpp"
#include "indexer.hpp"
#include "location.hpp"
#include "silo_output.hpp"
#include "node.hpp"
#include "tree.hpp"

#endif /* XTREE_HPP_ */
