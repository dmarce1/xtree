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
		this_range.first = 0.0;
		this_range.second = 1.0;
		this_sum = 0.0;
	}
	double func(double x) {
		return x * x;
	}
	void setter_a(xtree::location<1> requester, std::pair<double, double> range) {
		//	printf("calling setter_a\n");
		this_range = range;
		if (node_ptr->is_terminal()) {
			auto mid = (range.first + range.second) * 0.5;
			this_sum = (func(range.first) + 4.0 * func(mid) + func(range.second)) * (range.second - range.first) * (1.0 / 6.0);
		}
	}

	std::pair<double, double> getter_a(xtree::location<1> requester) {
		//	printf("calling getter_a\n");
		std::pair<double, double> range;
		if (requester.this_child_index() == 0) {
			range.first = this_range.first;
			range.second = (this_range.first + this_range.second) * 0.5;
		} else {
			range.first = (this_range.first + this_range.second) * 0.5;
			range.second = this_range.second;
		}
		return range;
	}

	void setter_d(xtree::location<1> requester, double sum) {
		//	printf("calling setter_d\n");
		this_sum += sum;
		if (node_ptr->get_level() == 0) {
			printf("%e\n", this_sum);
		}
	}

	bool getter_refine(xtree::location<1> self) {
		//	printf("getter_refine\n");
		if (self.get_level() < 10) {
			return true;
		} else {
			return false;
		}
	}

	double getter_d(xtree::location<1> requester) {
		//	printf("calling getter_d\n");
		return this_sum;
	}

	int exchange_test1() {
		return 1;
	}
	void exchange_test2(int) {
	}
};

int hpx_main() {
	hpx::id_type root_gid = (hpx::new_<xtree::tree_type>(hpx::find_here())).get();
	using op_refine = xtree::tree_type::operation<xtree::REBRANCH, bool, &IntegrationTree::getter_refine, nullptr>;
	using op_a = xtree::tree_type::operation<xtree::ASCEND, std::pair<double,double>, &IntegrationTree::getter_a, &IntegrationTree::setter_a>;
	using op_d = xtree::tree_type::operation<xtree::DESCEND, double, &IntegrationTree::getter_d, &IntegrationTree::setter_d>;
//	test.branch();
//	test.do_output();
	auto f1 = hpx::async<xtree::tree_type::action_execute_operators<std::tuple<op_refine>>>(root_gid);
	f1.get();
	for (int i = 1; i < 8; i++) {
		f1 = hpx::async<xtree::tree_type::action_execute_operators<std::tuple<op_refine>>>(root_gid);
		f1.get();
	}
	auto f2 = hpx::async<xtree::tree_type::action_execute_operators<std::tuple<op_a, op_d>>>(root_gid);
	f2.get();
	return hpx::finalize();
}
