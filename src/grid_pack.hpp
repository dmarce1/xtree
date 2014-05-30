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

template<typename GridBase, typename Tuple, int Iter = (std::tuple_size<Tuple>::value - 1)>
struct make_grid_pack {
	using grid_array_type = std::array<std::shared_ptr<GridBase>, std::tuple_size<Tuple>::value>;
	void operator()(grid_array_type& grids) const {
		make_grid_pack<GridBase, Tuple, Iter - 1> make_more;
		grids[Iter] = std::dynamic_pointer_cast<GridBase>(std::tuple_element<Iter, Tuple>::type::create());
		make_more(grids);
	}
};

template<typename GridBase, typename Tuple>
struct make_grid_pack<GridBase, Tuple, 0> {
	using grid_array_type = std::array<std::shared_ptr<GridBase>, std::tuple_size<Tuple>::value>;
	void operator()(grid_array_type& grids) const {
		grids[0] = std::dynamic_pointer_cast<GridBase>(std::tuple_element<0, Tuple>::type::create());
	}
};

template<typename ...Params>
class grid_pack {
	using grid_tuple = std::tuple<Params...>;
	using first_grid_type = typename std::tuple_element<0,grid_tuple>::type;
	using grid_base_type = grid_base<typename first_grid_type::type, first_grid_type::Ndim>;
	using grid_array_type = std::array<std::shared_ptr<grid_base_type>,sizeof...(Params)>;
	static constexpr int Ndim = first_grid_type::Ndim;
	grid_array_type grids;
	const node<grid_pack<Params...>, Ndim>& node_ref;
public:
	grid_pack(const node<grid_pack<Params...>, Ndim>& nref) :
			node_ref(nref) {
		make_grid_pack<grid_base_type, grid_tuple> make;
		make(grids);
	}
	template<int N>
	using grid_type = typename std::tuple_element<N, grid_tuple>::type;
	template<int N>
	using Dims = typename grid_type<N>::dims;
	template<int N>
	using T = typename grid_type<N>::type;
	template<int N>
	using descend_type = vector<typename grid_type<N>::type, grid_type<N>::Size / grid_type<N>::Nchild>;
	template<int N>
	std::shared_ptr<grid_type<N>> get_grid() {
		return std::static_pointer_cast < grid_type < N >> (grids[N]);
	}
	template<int N>
	descend_type<N> get_descend(const location<Ndim>& requester) {
		return get_grid<N>()->restrict_to_stream();
	}
	template<int N>
	void set_descend(const location<Ndim>& sender, descend_type<N> stream) {
		get_grid<N>()->restrict_from_stream(sender.this_child_index(), stream);
	}
	template<int N>
	bgrid<T<N>, Dims<N>, grid_type<N>::bw> get_decomp_boundary(const location<Ndim>& requester) {
		return bgrid<T<N>, Dims<N>, grid_type<N>::bw>(requester.relative_direction_to(this->node_ref.get_self(), get_grid<N>()));
	}
	template<int N>
	void set_decomp_boundary(const location<Ndim>& sender, bgrid<T<N>, Dims<N>, grid_type<N>::bw> bg) {
		get_grid<N>()->set_neighbor(std::make_shared<bgrid<T<N>, Dims<N>, grid_type<N>::bw>>(std::move(bg)),
				this->node_ref.get_self().relative_direction_to(sender));
	}
	template<int N>
	typename agrid<T<N>, Dims<N>, grid_type<N>::bw>::neighbors_type get_amr_boundary(const location<Ndim>& receiver) {
		typename agrid<T<N>, Dims<N>, grid_type<N>::bw>::neighbors_type v(pow_<3, Ndim>);
		int cnt = 0;
		for (dir_type<Ndim> dir; !dir.end(); dir++) {
			if (node_ref.get_boundary_type(dir) == AMR) {
				v[cnt].first = dir;
				v[cnt].second = agrid<T<N>, Dims<N>, grid_type<N>::bw>(receiver.this_child_index(), dir, get_grid<N>());
				cnt++;
			}
		}
		v.resize(cnt);
		return v;
	}
	template<int N>
	void set_amr_boundary(const location<Ndim>& sender, typename agrid<T<N>, Dims<N>, grid_type<N>::bw>::neighbors_type v) {
		for (int i = 0; i < v.size(); i++) {
			get_grid<N>()->set_neighbor(v[i].second, v[i].first);
		}
	}
};

}

#endif /* GRID_PACK_HPP_ */
