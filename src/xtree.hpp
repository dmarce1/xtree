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

#define XTREE_INSTANTIATE( ... )																							\
namespace xtree{																												\
typedef node<int_seq<__VA_ARGS__>> node_type;																								\
}																																\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::node_type>, node_type);		           	   	\
/**/

#endif /* XTREE_HPP_ */
