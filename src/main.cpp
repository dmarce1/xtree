/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include "xtree.hpp"

#define NDIM 3

struct mynode : xtree::node<mynode,NDIM>, hpx::components::simple_component_base<mynode> {
	typedef mynode type_holder;
    typedef xtree::node<mynode,NDIM> base_type_holder;

};

XTREE_INSTANTIATE(mynode, 3);

int hpx_main() {

	hpx::id_type tree_gid = (hpx::new_<xtree::tree_type>(hpx::find_here())).get();
	return hpx::finalize();
}
