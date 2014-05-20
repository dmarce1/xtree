/*
 * locality_server.hpp
 *
 *  Created on: May 20, 2014
 *      Author: dmarce1
 */

#ifndef LOCALITY_SERVER_HPP_
#define LOCALITY_SERVER_HPP_

#include "xtree.hpp"
#include "vector.hpp"

namespace xtree {
namespace server {
hpx::future<hpx::id_type> increment_load();
int get_load();
void decrement_load();
}
}

#endif /* LOCALITY_SERVER_HPP_ */
