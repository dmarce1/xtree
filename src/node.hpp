/*
 * node.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef NODE_HPP_
#define NODE_HPP_

namespace xtree {

template<typename, int, typename, op_type>
class ___setup_op_dataflow;

template<int Ndim>
class node_base {
public:
	node_base() {
	}
	virtual ~node_base() {
	}
	virtual const location<Ndim>& get_self() const = 0;
};

template<typename Member, int Ndim>
class node: public hpx::components::managed_component_base<node<Member, Ndim>>, public node_base<Ndim> {
public:
	static constexpr int Nchild = pow_<2, Ndim>::value;
	static constexpr int Nneighbor = pow_<3, Ndim>::value;
	static constexpr int Nniece = Nchild * Nneighbor;
	using dims_type = typename Member::dims_type;
	using base_type = hpx::components::managed_component_base<node<Member, Ndim>>;
	using child_array_type = std::array<hpx::id_type, Nchild>;
	using niece_array_type = std::array<std::array<hpx::id_type, Nchild>, Nneighbor>;
	using neighbor_notify_type = std::array<std::pair<int, hpx::id_type>, Nneighbor>;
	using neighbor_array_type = std::array<hpx::id_type, Nneighbor>;
private:
	static bool child_is_niece_of(child_index_type<Ndim> ci, dir_type<Ndim> dir) {
		for (int i = 0; i < Ndim; i++) {
			if (dir[i] == +1) {
				if (ci[i] == 0) {
					return false;
				}
			} else if (dir[i] == -1) {
				if (ci[i] == 1) {
					return false;
				}
			}
		}
		return true;
	}
private:
	bool is_leaf;
	child_array_type children;
	hpx::id_type parent;
	hpx::shared_future<void> flock, last_flock;
	location<Ndim> self;
	niece_array_type nieces;
	neighbor_array_type neighbors;
	std::shared_ptr<Member> data;
	tree<Member, Ndim>* local_tree;
private:
	hpx::id_type get_niece(child_index_type<Ndim> ci, dir_type<Ndim> dir) {
		child_index_type<Ndim> nci;
		dir_type<Ndim> ndi;
		for (int i = 0; i < Ndim; i++) {
			if (dir[i] == 0) {
				nci[i] = ci[i];
				ndi[i] = dir[i];
			} else if (dir[i] == +1) {
				ndi[i] = +(nci[i] ^ 0);
				nci[i] = nci[i] ^ 1;
			} else /*if( dir[i] == -1)*/{
				ndi[i] = -(nci[i] ^ 1);
				nci[i] = nci[i] ^ 1;
			}
		}
		return nieces[ndi][nci];
	}
public:
	node() {
		assert(false);
	}
	node(const location<Ndim>& _loc, hpx::id_type _parent_id, neighbor_array_type _neighbors, tree<Member, Ndim>* ltree) :
			local_tree(ltree) {
		data = std::make_shared < Member > (*this);
		self = _loc;
		flock = last_flock = hpx::make_ready_future();
		neighbors = _neighbors;
		is_leaf = true;
		for (int ni = 0; ni < Nneighbor; ni++) {
			for (int ci = 0; ci < Nchild; ci++) {
				nieces[ni][ci] = hpx::invalid_id;
			}
		}
		for (int ci = 0; ci < Nchild; ci++) {
			children[ci] = hpx::invalid_id;
		}
		parent = _parent_id;

	}

	virtual ~node() {
		if (!is_leaf) {
			debranch().get();
		}
		if (get_level() != 0) {
			local_tree->delete_node(this);
		}
	}

	template<typename Arc>
	void serialize(Arc& ar, const int) {
	}

	const location<Ndim>& get_self() const {
		return self;
	}

	bool is_terminal() const {
		return is_leaf;
	}

	hpx::future<void> branch() {
		if (!is_leaf) {
			return hpx::make_ready_future();
		}
		is_leaf = false;
		auto fut0 = hpx::make_ready_future();

		/* Allocate new nodes on localities */
		auto fut1 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](void) {
			std::vector<hpx::future<hpx::id_type>> futures(Nchild);
			for (child_index_type<Ndim> ci; !ci.end(); ci++) {
				neighbor_array_type pack;
				for( dir_type<Ndim> dir; !dir.end(); dir++ ) {
					pack[dir] = get_niece(ci,dir);
				}
				futures[ci] = local_tree->new_node(self.get_child(ci), base_type::get_gid(), pack);
			}
			return futures;
		}), std::move(fut0));

		/* Assign gids to children[], notify neighbors and children */
		auto fut2 = hpx::lcos::local::dataflow(hpx::util::unwrapped([this](std::vector<hpx::future<hpx::id_type>> vfut) {

			for (int ci = 0; ci < Nchild; ci++) {
				children[ci] = vfut[ci].get();
			}

			std::vector<hpx::future<void>> futures(Nneighbor);
			for( dir_type<Ndim> si; !si.end(); si++) {
				if( neighbors[si] != hpx::invalid_id ) {
					futures[si] = hpx::async<action_notify_branch>(neighbors[si], si, children);
				} else {
					futures[si] = hpx::make_ready_future();
				}
			}

			return futures;

		}), std::move(fut1));

		return when_all(fut2);
	}

	hpx::future<void> debranch() {
		is_leaf = true;
		std::vector<hpx::future<void>> nfutures(Nneighbor);
		std::vector<hpx::future<void>> cfutures(Nchild);
		for (child_index_type<Ndim> ci; !ci.end(); ci++) {
			if (children[ci] != hpx::invalid_id) {
				cfutures[ci] = hpx::async<action_debranch>(children[ci]);
			} else {
				cfutures[ci] = hpx::make_ready_future();
			}
		}
		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			if (neighbors[dir] != hpx::invalid_id) {
				nfutures[dir] = hpx::async<action_notify_debranch>(neighbors[dir], dir);
			} else {
				nfutures[dir] = hpx::make_ready_future();
			}
		}
		auto fut1 = when_all(cfutures).then(hpx::util::unwrapped([this](std::vector<hpx::future<void>>) {
			for (child_index_type<Ndim> ci; !ci.end(); ci++) {
				children[ci] = hpx::invalid_id;
			}
		}));
		fut1.get();
		return when_all(nfutures);
	}
	int get_level() const {
		return self.get_level();
	}

	hpx::future<void> notify_branch(dir_type<Ndim> dir, const child_array_type& nephews) {
		hpx::future<void> ret_fut;
		dir.flip();
		nieces[dir] = nephews;
		if (!is_leaf) {
			std::vector<hpx::future<void>> futs(Nchild / 2);
			int index = 0;
			for (child_index_type<Ndim> ci; !ci.end(); ci++) {
				if (child_is_niece_of(ci, dir)) {
					futs[index++] = hpx::async<action_notify_of_neighbor>(children[ci], dir, get_niece(ci, dir));
				}
			}
			futs.resize(index);
			ret_fut = when_all(std::move(futs));
		} else {
			ret_fut = hpx::make_ready_future();
		}
		return ret_fut;
	}

	hpx::future<void> notify_debranch(dir_type<Ndim> dir) {
		hpx::future<void> ret_fut;
		dir.flip();
		std::fill(nieces[dir].begin(), nieces[dir].end(), hpx::invalid_id);
		if (!is_leaf) {
			std::vector<hpx::future<void>> futs(Nchild / 2);
			int index = 0;
			for (child_index_type<Ndim> ci; !ci.end(); ci++) {
				if (child_is_niece_of(ci, dir)) {
					futs[index++] = hpx::async<action_notify_of_neighbor>(children[ci], dir, hpx::invalid_id);
				}
			}
			futs.resize(index);
			ret_fut = when_all(futs);
		} else {
			ret_fut = hpx::make_ready_future();
		}
		return ret_fut;
	}

	void notify_of_neighbor(dir_type<Ndim> dir, hpx::id_type id) {
		neighbors[dir] = id;
	}

	node<Member, Ndim>* get_this() {
		return this;
	}

	bound_type get_boundary_type(const dir_type<Ndim>& dir) const {
		if (neighbors[dir] != hpx::invalid_id) {
			return DECOMP;
		} else {
			for (int di = 0; di < Ndim; di++) {
				if (((dir[di] == -1) && (self.get_location(di) == 0)) || ((dir[di] == +1) && (self.get_location(di) == ((1 << self.get_level()) - 1)))) {
					return PHYS;
				}
			}
			return AMR;
		}
	}

	bool has_amr_boundary() const {
		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			if (get_boundary_type(dir) == AMR) {
				return true;
			}
		}
		return false;
	}

	template<typename, int, typename, op_type>
	friend class ___setup_op_dataflow;

	template<typename Op>
	using setup_op_dataflow = ___setup_op_dataflow<Member,Ndim,Op,Op::op>;

	hpx::shared_future<void> get_last_future() {
		return flock;
	}

	template<typename Op>
	hpx::future<typename Op::type> get_op_data(location<Ndim> requester) {
		return ((Op::op == EXCHANGE) ? last_flock : flock).then(hpx::util::unwrapped([this, requester]() {
			Op op;
			return op.get(data,requester);
		}));
	}
	template<typename Op>
	using action_get_op_data = typename hpx::actions::make_action<decltype(&node<Member,Ndim>::get_op_data<Op>), &node<Member,Ndim>::get_op_data<Op>>;

	XTREE_MAKE_ACTION(action_debranch, node::debranch);

	XTREE_MAKE_ACTION(action_get_this, node::get_this);

	XTREE_MAKE_ACTION(action_notify_branch, node::notify_branch);

	XTREE_MAKE_ACTION(action_notify_debranch, node::notify_debranch);

	XTREE_MAKE_ACTION(action_notify_of_neighbor, node::notify_of_neighbor);

}
;

template<typename Member, int Ndim, typename Op, op_type Optype>
struct ___setup_op_dataflow {
	void operator()(node<Member, Ndim>& node_ref) const {
		using action = typename node<Member, Ndim>::template action_get_op_data<Op>;
		hpx::shared_future<void> future;
		std::vector<hpx::shared_future<void>> futures;
		switch (Optype) {
		case AMR_ASCEND:
			if ((node_ref.parent != hpx::invalid_id) && node_ref.has_amr_boundary()) {
				auto fut = hpx::async<action>(node_ref.parent, node_ref.self);
				future = fut.then(hpx::util::unwrapped([&node_ref](typename Op::type arg) {
					Op op;
					op.set( node_ref.data, node_ref.self.get_parent(), arg);
				}));
			} else {
				future = hpx::make_ready_future();
			}
			break;
		case ASCEND:
			if (node_ref.parent != hpx::invalid_id) {
				auto fut = hpx::async<action>(node_ref.parent, node_ref.self);
				future = fut.then(hpx::util::unwrapped([&node_ref](typename Op::type arg) {
					Op op;
					op.set( node_ref.data, node_ref.self.get_parent(), arg);
				}));
			} else {
				future = hpx::make_ready_future();
			}
			break;
		case DESCEND:
			if (!node_ref.is_leaf) {
				futures.resize(node_ref.Nchild);
				for (child_index_type<Ndim> ci; !ci.end(); ci++) {
					auto fut = hpx::async<action>(node_ref.children[ci], node_ref.self);
					futures[ci] = fut.then(hpx::util::unwrapped([&](typename Op::type arg) {
						Op op;
						op.set( node_ref.data, node_ref.self.get_child(ci),arg);
					}));
				}
			} else {
				futures.resize(1);
				futures[0] = hpx::make_ready_future();
			}
			future = when_all(futures).share();
			break;
		case EXCHANGE:
			futures.resize(node_ref.Nneighbor);
			for (dir_type<Ndim> dir; !dir.end(); dir++) {
				if (node_ref.neighbors[dir] != hpx::invalid_id) {
					auto fut = hpx::async<action>(node_ref.children[dir], node_ref.self);
					futures[dir] = fut.then(hpx::util::unwrapped([&node_ref,dir](typename Op::type arg) {
						Op op;
						op.set( node_ref.data, node_ref.self.get_neighbor(dir),arg);
					}));
				} else {
					futures[dir] = hpx::make_ready_future();
				}
			}
			future = when_all(futures).share();
			break;
		}
		node_ref.last_flock = node_ref.flock;
		node_ref.flock = when_all(future, node_ref.last_flock).share();
	}
};

template<typename Member, int Ndim, typename Op>
struct ___setup_op_dataflow<Member, Ndim, Op, op_type::REBRANCH> {
	void operator()(node<Member, Ndim>& node_ref) const {
		using action = typename node<Member, Ndim>::template action_get_op_data<Op>;
		hpx::shared_future<void> future;
		if (node_ref.is_leaf) {
			Op op;
			if (op.get(node_ref.data, node_ref.self)) {
				future = node_ref.branch().share();
			} else {
				node_ref.debranch();
				future = hpx::make_ready_future().share();
			}
		} else {
			future = hpx::make_ready_future().share();
		}
		node_ref.last_flock = node_ref.flock;
		node_ref.flock = when_all(future, node_ref.last_flock).share();
	}
};

template<typename Member, int Ndim, typename Op>
struct ___setup_op_dataflow<Member, Ndim, Op, op_type::LOCAL> {
	void operator()(node<Member, Ndim>& node_ref) const {
		hpx::shared_future<void> future = (node_ref.flock).then(hpx::util::unwrapped([&]() {
			Op op;
			op.get(node_ref.data, node_ref.self);
		}));
		node_ref.last_flock = node_ref.flock;
		node_ref.flock = when_all(future, node_ref.last_flock).share();
	}

};

}

#endif /* NODE_HPP_ */
