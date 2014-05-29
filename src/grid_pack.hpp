/*
 * grid_pack.hpp
 *
 *  Created on: May 29, 2014
 *      Author: dmarce1
 */

#ifndef GRID_PACK_HPP_
#define GRID_PACK_HPP_

#include "fwd.hpp"

namespace xtree {

template<typename Grid, typename ...Params>
class grid_pack: public grid_pack<Params...> {
public:
	using T = typename Grid::type;
	using restrict_type = vector<T, Grid::Size / Grid::Nchild>;
	using Dims = typename Grid::dims_type;
	static constexpr int Bw = Grid::bw;
	static constexpr int Ndim = Dims::dim();
	using base_type = grid_pack<Params...>;
private:
	std::shared_ptr<Grid> this_grid;
public:
	template<int N>
	restrict_type get_descend_restrict(const location<Ndim>& requester) {
		if (N == 0) {
			return this_grid->restrict_to_stream();
		} else {
			return base_type::get_descend_restrict(requester);
		}
	}
	template<int N>
	void set_descend_restrict(const location<Ndim>& sender, restrict_type stream) {
		if (N == 0) {
			this_grid->restrict_from_stream(sender.this_child_index(), stream);
		} else {
			base_type::restrict_from_stream(sender, stream);
		}
	}
	template<int N>
	bgrid<T, Dims, Bw> get_decomp_boundary(const location<Ndim>& requester) {
		if (N == 0) {
			return bgrid<T, Dims, Bw>(requester.relative_direction_to(this->node_ptr->get_self()));
		} else {
			return base_type::get_decomp_boundary<N - 1>(requester);
		}
	}
	template<int N>
	void set_decomp_boundary(const location<Ndim>& sender, bgrid<T, Dims, Bw> bg) {
		if (N == 0) {
			this_grid->set_neighbor(std::make_shared<bgrid<T, Dims, Bw>>(std::move(bg)), this->node_ptr->get_self().relative_direction_to(sender));
		} else {
			base_type::set_decomp_boundary<N - 1>(sender, bg);
		}
	}

	grid_pack(const node_base<Ndim>* nodeptr) :
			base_type(nodeptr) {
		this_grid = Grid::create();
	}
	virtual ~grid_pack() {
	}

};

template<typename Grid>
class grid_pack<Grid> {
public:
	using T = typename Grid::type;
	using restrict_type = vector<T, Grid::Size / Grid::Nchild>;
	using Dims = typename Grid::dims_type;
	static constexpr int Bw = Grid::bw;
	static constexpr int Ndim = Dims::dim();
private:
	std::shared_ptr<Grid> this_grid;
protected:
	const node_base<Ndim>* node_ptr;
public:
	template<int N>
	restrict_type get_descend_restrict(const location<Ndim>& requester) {
		return this_grid->restrict_to_stream();
	}
	template<int N>
	void set_descend_restrict(const location<Ndim>& sender, restrict_type stream) {
		this_grid->restrict_from_stream(sender.this_child_index(), stream);
	}
	template<int N>
	bgrid<T, Dims, Bw> get_decomp_boundary(const location<Ndim>& requester) {
		return bgrid<T, Dims, Bw>(requester.relative_direction_to(node_ptr->get_self()));
	}
	template<int N>
	void set_decomp_boundary(const location<Ndim>& sender, bgrid<T, Dims, Bw> bg) {
		this_grid->set_neighbor(std::make_shared<bgrid<T, Dims, Bw>>(std::move(bg)), node_ptr->get_self().relative_direction_to(sender));
	}
	grid_pack(const node_base<Ndim>* nodeptr) :
			node_ptr(nodeptr) {
		this_grid = Grid::create();
	}
	virtual ~grid_pack() {
	}

}
;
}

#endif /* GRID_PACK_HPP_ */
