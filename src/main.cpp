/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include <mpi.h>
#include "xtree.hpp"
#include "./fmmx_node.hpp"
//#include <fenv.h>

using namespace xtree;
using namespace fmmx;




using fmmx_node_type = fmmx_node<3,8,5>;
XTREE_INSTANTIATE(fmmx_node_type, 3);

using rg_func = fmmx_node_type::regrid_function<&fmmx_node_type::regrid_test>;
using init_func = fmmx_node_type::local_function<&fmmx_node_type::init_grid>;
using dfunc = fmmx_node_type::descend_function<fmmx_node_type::descend_type, &fmmx_node_type::descend>;
using eg_func = fmmx_node_type::exchange_get_function<fmmx_node_type::exchange_type, &fmmx_node_type::exchange_get>;
using es_func = fmmx_node_type::exchange_set_function<fmmx_node_type::exchange_type, &fmmx_node_type::exchange_set>;
using afunc = fmmx_node_type::ascend_function<fmmx_node_type::ascend_type, &fmmx_node_type::ascend>;

typedef fmmx_node_type::action_regrid<rg_func> action_regrid;
typedef fmmx_node_type::action_local<init_func> action_local;
typedef fmmx_node_type::action_descend<dfunc> action_descend_func;
typedef fmmx_node_type::action_ascend<afunc> action_ascend_func;
typedef fmmx_node_type::action_exchange_get<eg_func,es_func> action_get_func;
typedef fmmx_node_type::action_exchange_set<es_func> action_set_func;


HPX_REGISTER_ACTION( action_regrid );
HPX_REGISTER_ACTION( action_local );
HPX_REGISTER_ACTION( action_descend_func );
HPX_REGISTER_ACTION( action_ascend_func );
HPX_REGISTER_ACTION( action_get_func );
HPX_REGISTER_ACTION( action_set_func);

void test() {
}


void execute() {
	printf( "Initializing program...\n");
	hpx::id_type tree_gid = (hpx::new_<tree_type>(hpx::find_here())).get();
		auto fut0 = hpx::async<tree_type::action_get_this>(tree_gid);
	tree_type* tree_ptr = fut0.get();
	tree_ptr->place_root();
	auto fut1 = tree_ptr->get_root();
	fmmx_node_type* root_node = fut1.get();
	std::vector<fmmx_node_type::operation_type> refine_ops(1);
	std::vector<fmmx_node_type::operation_type> init_ops(1);
	refine_ops[0] = fmmx_node_type::make_regrid_operation<rg_func>();
	init_ops[0] = fmmx_node_type::make_local_operation<init_func>();
	printf( "Refining...\n");
	root_node->regrid<rg_func>(0);
	printf( "Refining...\n");
	root_node->regrid<rg_func>(0);
	printf( "Refining...\n");
	root_node->regrid<rg_func>(0);
	printf( "Refining...\n");
	root_node->regrid<rg_func>(0);
	printf( "Initializing grid...\n");
	root_node->execute_operations(init_ops);

		std::vector<fmmx_node_type::operation_type> ops(3);
	auto dop = fmmx_node_type::make_descend_operation<dfunc>();
	auto eop = fmmx_node_type::make_exchange_operation<eg_func, es_func>();
	auto aop = fmmx_node_type::make_ascend_operation<afunc>();
	ops[0] = std::move(dop);
	ops[1] = std::move(eop);
	ops[2] = std::move(aop);
	double tstart = MPI_Wtime();
	printf( "Starting \n");
	root_node->execute_operations(ops);

	printf( "Done in %e seconds\n", MPI_Wtime() - tstart );
	printf("Output\n");
	tree_ptr->output();
		printf( "Output done\n");
	tree_ptr->delete_node(root_node);
}

int hpx_main() {
	printf( "HPX MAIN\n");
#ifndef NDEBUG
//	feenableexcept(FE_DIVBYZERO);
//	feenableexcept(FE_INVALID);
//	feenableexcept(FE_OVERFLOW);
#endif
	execute();
	return hpx::finalize();

}
