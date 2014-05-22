/*
 * xtree.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef XTREE_HPP_
#define XTREE_HPP_

#include "int_seq.hpp"
#include "node.hpp"

#define XTREE_INSTANTIATE( DERIVED_CLASS, ... )																							\
namespace xtree{																												\
typedef node<DERIVED_CLASS, xtree::int_seq<__VA_ARGS__>> base_node_type;   \
}																																\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::base_node_type>, base_node_type);		           	   	\
/**/

#endif /* XTREE_HPP_ */
