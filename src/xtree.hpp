/*
 * xtree.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef XTREE_HPP_
#define XTREE_HPP_

#include <hpx/hpx_init.hpp>
#include <hpx/include/components.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/lcos/local/counting_semaphore.hpp>
#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/util/unwrapped.hpp>

#include <boost/mpl/int.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>

#include "indexer.hpp"
#include "int_seq.hpp"
#include "locality_server.hpp"
#include "node.hpp"
#include "pow_.hpp"
#include "vector.hpp"

#define XTREE_INSTANTIATE( __NDIM__ )																							\
namespace xtree{																												\
typedef node<__NDIM__> node_type;																								\
}																																\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::node_type>, node_type);		           	   	\
/**/

#endif /* XTREE_HPP_ */
