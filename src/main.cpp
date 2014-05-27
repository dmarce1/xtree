/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include "xtree.hpp"

class IntegrationTree;

XTREE_INSTANTIATE(IntegrationTree, 1);

class IntegrationTree {
private:
	std::pair<double, double> this_range;
	double this_sum;
	const xtree::node<IntegrationTree, 1>* node_ptr;
public:
	IntegrationTree(const xtree::node<IntegrationTree, 1>* ptr) {
		node_ptr = ptr;
	}
	double func(double x) {
		return x * x;
	}
	int descend_test(std::vector<int> inputs) {
		return 0;
	}
	std::vector<int> ascend_test(int input) {
		return std::vector<int>(0);
	}

	void setter_a(xtree::location<1> requester, std::pair<double, double> range) {
		this_sum = 0.0;
		this_range = range;
		if (node_ptr->is_terminal()) {
			this_sum = (func(range.first) + func(range.second)) * (range.second - range.first);
		}
	}

	std::pair<double, double> getter_a(xtree::location<1> requester) {
		std::pair<double, double> range;
		if (requester.this_child_index() == 0) {
			range.first = this_range.first;
			range.second = (this_range.first + this_range.second) * 0.5;
		} else {
			range.first = (this_range.first + this_range.second) * 0.5;
			range.second = this_range.first;
		}
		return range;
	}

	void setter_d(xtree::location<1> requester, double sum) {
		this_sum += sum;
	}

	double getter_d(xtree::location<1> requester) {
		return this_sum;
	}

	int exchange_test1() {
		return 1;
	}
	void exchange_test2(int) {
	}
};

int hpx_main() {
	xtree::tree_type test;
	using op_a = xtree::tree_type::operation<xtree::ASCEND, std::pair<double,double>, &IntegrationTree::getter_a, &IntegrationTree::setter_a>;
	using op_d = xtree::tree_type::operation<xtree::DESCEND, double, &IntegrationTree::getter_d, &IntegrationTree::setter_d>;
//	test.branch();
//	test.do_output();
	test.execute_operators<std::tuple<op_a, op_d>>();
	return hpx::finalize();
}
