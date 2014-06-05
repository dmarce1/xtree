/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include "xtree.hpp"
/*
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
 };*/

//using namespace xtree;
class mynode;
typedef xtree::grid<double, xtree::int_seq<8, 8, 8>, 1> grid_type1;
//typedef xtree::grid<int, xtree::int_seq<8, 8, 8>, 2> grid_type2;
typedef xtree::grid_pack<mynode, grid_type1> grid_pack_type;

class mynode: public grid_pack_type {
public:
	mynode() = delete;
	mynode(const xtree::node<mynode, 3>& nref) :
			grid_pack_type(nref) {
	}
	struct refine_check {
		bool operator()(const std::shared_ptr<mynode>& _this, const xtree::location<Ndim>&) const {
			if (_this->node_ref.get_level() < 2) {
				return true;
			} else {
				return false;
			}
		}
	};
	template<typename T>
	struct null_set {
		void operator()(const std::shared_ptr<mynode>&, const xtree::location<Ndim>&, const T&) const {
		}
	};
	virtual ~mynode() {
	}
};
XTREE_INSTANTIATE(mynode, 3);

template<int N>
using descend_operation = std::tuple<
xtree::tree_type::operation<xtree::DESCEND,
grid_pack_type::descend_type<N>,
mynode::get_descend<N>,
mynode::set_descend<N>
>>;
using refine_operation = std::tuple<xtree::tree_type::operation<xtree::REBRANCH,bool,mynode::refine_check,mynode::null_set<bool>>>;

int hpx_main() {
	hpx::id_type tree_gid = (hpx::new_<xtree::tree_type>(hpx::find_here())).get();
	auto f0 = hpx::async<xtree::tree_type::action_place_root>(tree_gid);
	f0.get();
	auto f1 = hpx::async<xtree::tree_type::action_execute_operators<refine_operation>>(tree_gid);
	f1.get();
	auto f2 = hpx::async<xtree::tree_type::action_output>(tree_gid);
	f2.get();
	return hpx::finalize();
}
