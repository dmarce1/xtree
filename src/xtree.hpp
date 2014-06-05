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

#define XTREE_MAKE_ACTION( a, b ) 																					\
using a = typename hpx::actions::make_action<decltype(&b),&b>::type;	 											\
/**/

#define XTREE_INSTANTIATE( MEMBER_CLASS, ... )																		\
namespace xtree {																									\
	typedef silo_output<__VA_ARGS__, 0> silo_output_type;   														\
	typedef node<MEMBER_CLASS, __VA_ARGS__> node_type;   														\
	typedef tree<MEMBER_CLASS, __VA_ARGS__> tree_type;   															\
}																													\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::silo_output_type>, silo_output_type);	\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::tree_type>, tree_type);	\
HPX_REGISTER_MINIMAL_COMPONENT_FACTORY(hpx::components::managed_component<xtree::node_type>, node_type);	\
/**/


#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/lcos/when_all.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/util/unwrapped.hpp>

#include <boost/mpl/int.hpp>
#include <boost/serialization/vector.hpp>

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
void serialize(Archive & ar, std::array<T,N> & a, const unsigned int version)
{
  ar & boost::serialization::make_array(a.data(), a.size());
}

} // namespace serialization
} // namespace boost

#include "load_balancer.hpp"
#include "util.hpp"
#include "vector.hpp"
#include "indexer.hpp"
#include "location.hpp"
#include "grid_index.hpp"
#include "grid_base.hpp"
#include "grid.hpp"
#include "agrid.hpp"
#include "bgrid.hpp"
#include "xgrid.hpp"
#include "grid_pack.hpp"
#include "silo_output.hpp"
#include "node.hpp"
#include "tree.hpp"
#include "state.hpp"


#endif /* XTREE_HPP_ */
