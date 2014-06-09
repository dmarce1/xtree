/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include "xtree.hpp"

#include <boost/multi_array.hpp>

namespace fmmx {

template<int Ndim>
struct center_of_mass {
	double mass;
	std::array<double, Ndim> position;
	template<typename Arc>
	void serialize(Arc& ar, const int v) {
		ar & mass;
		serialize(ar, position, v);
	}
};

template<int Ndim, int Nx>
class gnode: public xtree::node<gnode<Ndim, Nx>, Ndim>, public hpx::components::managed_component_base<gnode<Ndim, Nx>> {
public:
	using type_holder = gnode;
	using base_type_holder = xtree::node<gnode, Ndim>;
	using wrapped_type = gnode;
	using com_array_type = boost::multi_array<center_of_mass<Ndim>, Ndim>;
private:
	com_array_type com_array;
public:
	com_array_type com_descend(const std::vector<com_array_type>& child_arrays) {
		com_array_type this_array;
		return this_array;
	}
	std::vector<double> L_ascend(const double& ds) {
		return std::vector<double>();
	}

};

}

using gnode_type = fmmx::gnode<3,8>;
XTREE_INSTANTIATE(gnode_type, 3);

int hpx_main() {
	hpx::id_type root_gid = (hpx::new_<xtree::node_type>(hpx::find_here())).get();
	xtree::tree_type* tree_ptr = hpx::async<xtree::tree_type::action_get_this>(root_gid).get();
//xtree::node_type::operation_type op1 =
	typedef xtree::node_type::descend_function<gnode_type::com_array_type, &gnode_type::com_descend> descend_func;
	typedef xtree::node_type::ascend_function<double, &gnode_type::L_ascend> ascend_func;
	gnode_type::make_descend_operation<gnode_type::com_array_type, descend_func>();
	gnode_type::make_ascend_operation<double, ascend_func>();
//	node_type* root_node = tree_ptr->get_root().get();
	return hpx::finalize();
}
