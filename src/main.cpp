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

	hpx::id_type tree_gid = (hpx::new_<tree_type>(hpx::find_here())).get();
	auto fut0 = hpx::async<tree_type::action_get_this>(tree_gid);
	tree_type* tree_ptr = fut0.get();
	tree_ptr->place_root();
	auto fut1 = tree_ptr->get_root();
	fmmx_node_type* root_node = fut1.get();
	using dfunc = fmmx_node_type::descend_function<fmmx_node_type::descend_type, &fmmx_node_type::descend>;
	using rg_func = fmmx_node_type::regrid_function<&fmmx_node_type::regrid_test>;
	using eg_func = fmmx_node_type::exchange_get_function<fmmx_node_type::exchange_type, &fmmx_node_type::exchange_get>;
	using es_func = fmmx_node_type::exchange_set_function<fmmx_node_type::exchange_type, &fmmx_node_type::exchange_set>;
	using afunc = fmmx_node_type::ascend_function<fmmx_node_type::ascend_type, &fmmx_node_type::ascend>;

	std::vector<fmmx_node_type::operation_type> init_ops(1);
	init_ops[0] = fmmx_node_type::make_regrid_operation<rg_func>();
	root_node->execute_operations(init_ops);
	tree_ptr->output();
	/*
	 std::vector<fmmx_node_type::operation_type> ops(3);
	 auto dop = fmmx_node_type::make_descend_operation<dfunc>();
	 auto eop = fmmx_node_type::make_exchange_operation<eg_func, es_func>();
	 auto aop = fmmx_node_type::make_ascend_operation<afunc>();
	 ops[0] = dop;
	 ops[1] = eop;
	 ops[2] = aop;
	 root_node->execute_operations(ops);*/
	return hpx::finalize();
}
