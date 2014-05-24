/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include "xtree.hpp"

class Tree;

XTREE_INSTANTIATE(Tree, 8, 8, 8);

class Tree: public xtree::base_node_type {
public:
	int descend_test(std::vector<int> inputs) {
		return 0;
	}
	std::vector<int> ascend_test(int input) {
		return std::vector<int>(0);
	}

	void setter(xtree::location<Ndim> requester, int) {

	}

	int getter(xtree::location<Ndim> requester) {
		return 1;
	}

	int exchange_test1() {
		return 1;
	}
	void exchange_test2(int) {
	}
};

int hpx_main() {
	xtree::tree_type test;
	test.begin_execute<xtree::base_node_type::ASCEND, int, &Tree::getter, &Tree::setter>();
	return hpx::finalize();
}
