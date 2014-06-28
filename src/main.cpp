/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include "xtree.hpp"
#include "./fmmx_node.hpp"

using namespace xtree;
using namespace fmmx;

using fmmx_node_type = fmmx_node<3,8,5>;
XTREE_INSTANTIATE(fmmx_node_type, 3);

void test() {
}

int hpx_main() {

	hpx::id_type root_gid = (hpx::new_<node_type>(hpx::find_here())).get();
	tree_type* tree_ptr = hpx::async<tree_type::action_get_this>(root_gid).get();
	fmmx_node_type* root_node = tree_ptr->get_root().get();
	using dfunc = fmmx_node_type::descend_function<fmmx_node_type::descend_type, &fmmx_node_type::descend>;
	using eg_func = fmmx_node_type::exchange_get_function<fmmx_node_type::exchange_type, &fmmx_node_type::exchange_get>;
	using es_func = fmmx_node_type::exchange_set_function<fmmx_node_type::exchange_type, &fmmx_node_type::exchange_set>;
	using afunc = fmmx_node_type::ascend_function<fmmx_node_type::ascend_type, &fmmx_node_type::ascend>;
	std::vector<fmmx_node_type::operation_type> ops(3);
	auto dop = fmmx_node_type::make_descend_operation<dfunc>();
	auto eop = fmmx_node_type::make_exchange_operation<eg_func, es_func>();
	auto aop = fmmx_node_type::make_ascend_operation<afunc>();
	ops[0] = dop;
	ops[1] = eop;
	ops[2] = aop;
	root_node->execute_operations(ops);
	return hpx::finalize();
}
