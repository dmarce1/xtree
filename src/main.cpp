/*
 * main.cpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#include "xtree.hpp"

class mynode;
typedef xtree::grid<xtree::state<3>, xtree::int_seq<8, 8, 8>, 1> grid_type1;
typedef xtree::grid_pack<mynode, grid_type1> grid_pack_type;

class mynode: public grid_pack_type {
public:
	using dims_type = grid_type1::dims_type;
	static constexpr int Ndim = dims_type::dim();
private:
	std::array<double, Ndim> dx;
	grid_type1& Uc;
	double this_dtmin;
	double t;
public:
	mynode() = delete;
	mynode(const xtree::node<mynode, 3>& nref) :
			grid_pack_type(nref), Uc(*std::static_pointer_cast < grid_type1 > (grid_pack_type::grids[0])) {
		const double dx0 = nref.get_self().get_dx();
		for (int di = 0; di < Ndim; di++) {
			dx[di] = dx0 / dims_type::get(di);
		}
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
	struct null_set {
		void operator()(const std::shared_ptr<mynode>&, const xtree::location<Ndim>&) const {
		}
	};
	virtual ~mynode() {
	}
	struct local_operation {
		void operator()(const std::shared_ptr<mynode>&, const xtree::location<Ndim>& self) const {
			printf("Local operation\n");
		}
	};
	struct initialize_op {
		void operator()(const std::shared_ptr<mynode>& _this, const xtree::location<Ndim>& self) const {
			_this->this_dtmin = std::numeric_limits<double>::max();
		}
	};

	struct initialize_at_t0_op {
		void operator()(const std::shared_ptr<mynode>& _this, const xtree::location<Ndim>&) const {
			_this->initialize_at_t0();
		}
	};

	struct cfl_descend_get_op {
		double operator()(const std::shared_ptr<mynode>& _this, const xtree::location<Ndim>& self) const {
			_this->this_dtmin = std::min(0.4 / _this->get_maxdtinv(), _this->this_dtmin);
			printf("descend:get %e\n", _this->this_dtmin);
			return _this->this_dtmin;
		}
	};
	struct cfl_descend_set_op {
		void operator()(const std::shared_ptr<mynode>& _this, const xtree::location<Ndim>&, double mindt) const {
			_this->this_dtmin = std::min(_this->this_dtmin, mindt);
			printf("descend:set %e\n", _this->this_dtmin);
		}
	};

	struct cfl_ascend_get_op {
		double operator()(const std::shared_ptr<mynode>& _this, const xtree::location<Ndim>&) const {
			printf("ascend:get %e\n", _this->this_dtmin);
			return _this->this_dtmin;
		}
	};
	struct cfl_ascend_set_op {
		void operator()(const std::shared_ptr<mynode>& _this, const xtree::location<Ndim>&, double mindt) const {
			printf("ascend:set %e\n", _this->this_dtmin);
			_this->this_dtmin = std::min(_this->this_dtmin, mindt);
		}
	};

	friend class initialize_op;
	friend class initialize_at_t0_op;
	friend class cfl_descend_get_op;
	friend class cfl_descend_set_op;
	friend class cfl_ascend_get_op;
	friend class cfl_ascend_set_op;

	double get_maxdtinv() const {
		double max_dtinv = 0.0;
		for (xtree::grid_index<Ndim> gi(dims_type::to_vector()); !gi.end(); gi++) {
			for (int di = 0; di < Ndim; di++) {
				max_dtinv = std::max(max_dtinv, Uc[gi].get_signal_speed(di) / dx[di]);
			}
		}
		return max_dtinv;
	}

	std::array<double, Ndim> get_position(xtree::grid_index<Ndim>& gi) const {
		return node_ref.get_self().get_position() + dx * (gi.to_double() + 0.5);
	}

	double get_position(const xtree::grid_index<Ndim>& gi, int d) const {
		return node_ref.get_self().get_position()[d] + dx[d] * (double(gi[d]) + 0.5);
	}

	void initialize_at_t0() {
		printf("Computing CFL\n");
		t = 0.0;
		for (xtree::grid_index<Ndim> gi(dims_type::to_vector()); !gi.end(); gi++) {
			for (int di = 0; di < Ndim; di++) {
				Uc[gi].set_momentum(di, 0.0);
			}
			if (get_position(gi, 0) > 0.5) {
				Uc[gi].set_density(1.0);
				Uc[gi].set_energy(2.5);
			} else {
				Uc[gi].set_density(1.0e-1);
				Uc[gi].set_energy(1.25e-1);
			}
		}
	}

}
;
XTREE_INSTANTIATE(mynode, 3);

using init_t0_operation = xtree::tree_type::operation<xtree::LOCAL , void, mynode::initialize_at_t0_op>;
using init_operation = xtree::tree_type::operation<xtree::LOCAL , void, mynode::initialize_op>;
using cfl_descend_operation = xtree::tree_type::operation<xtree::DESCEND, double, mynode::cfl_descend_get_op, mynode::cfl_descend_set_op>;
using cfl_ascend_operation = xtree::tree_type::operation<xtree::ASCEND , double, mynode::cfl_ascend_get_op, mynode::cfl_ascend_set_op>;

using startup_operation = std::tuple<init_t0_operation >;
using hydro_step_operation = std::tuple<init_operation, cfl_descend_operation, cfl_ascend_operation >;
using refine_operation = std::tuple<xtree::tree_type::operation<xtree::REBRANCH,bool,mynode::refine_check>>;

int hpx_main() {

	hpx::id_type tree_gid = (hpx::new_<xtree::tree_type>(hpx::find_here())).get();
	auto f0 = hpx::async<xtree::tree_type::action_place_root>(tree_gid);
	f0.get();
	auto f1 = hpx::async<xtree::tree_type::action_execute_operators<refine_operation>>(tree_gid);
	f1.get();
	auto f4 = hpx::async<xtree::tree_type::action_execute_operators<startup_operation>>(tree_gid);
	f4.get();
	auto f2 = hpx::async<xtree::tree_type::action_output>(tree_gid);
	f2.get();
	auto f3 = hpx::async<xtree::tree_type::action_execute_operators<hydro_step_operation>>(tree_gid);
	f3.get();
	return hpx::finalize();
}
