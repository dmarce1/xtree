/*
 * array.hpp
 *
 *  Created on: May 27, 2014
 *      Author: dmarce1
 */

#ifndef ARRAY_HPP_
#define ARRAY_HPP_

#include "fwd.hpp"

namespace xtree {

template<typename T, typename Dims, int Bw = 1>
class grid: public grid_base<T, Dims::dim()> {
	friend class bgrid<T, Dims, Bw> ;
public:
	static constexpr int Ndim = Dims::dim();
	using base_type = grid_base<T, Ndim>;
	using index_type = typename grid_base<T, Ndim>::index_type;
	using dir_type = indexer<Ndim, 3, -1>;
private:
	static constexpr int Size = Dims::size();
	std::array<T, Size> data;
	std::array<const base_type*, pow_<3, Ndim>::value> grid_selector;
private:
	T& get(const index_type& i) {
		return data[base_type::template vector_to_index<Dims>(i)];
	}
	const T& get(const index_type& i) const {
		return data[base_type::template vector_to_index<Dims>(i)];
	}
public:
	grid() {
		for (dir_type i; !i.end(); i++) {
			grid_selector[i] = nullptr;
		}
		dir_type i;
		for (int di = 0; di < Ndim; di++) {
			i[di] = 0;
		}
		grid_selector[i] = this;
	}
	T& operator[](const index_type& i) {
		return get(i);
	}
	const T& operator[](const index_type& global_index) const {
		dir_type grid_index;
		index_type local_index;
		for (int i = 0; i < Ndim; i++) {
			const div_t tmp1 = div(global_index[i] + Dims::get(i), Dims::get(i));
			grid_index[i] = tmp1.quot - 1;
			local_index[i] = tmp1.rem;
		}
		return grid_selector[grid_index]->get(local_index);
	}
	void set_neighbor(const base_type* n, dir_type dir) {
		grid_selector[dir] = n;
	}
	virtual ~grid() {
	}
};

}

#endif /* ARRAY_HPP_ */
