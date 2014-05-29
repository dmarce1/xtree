/*
 * bgrid.hpp
 *
 *  Created on: May 29, 2014
 *      Author: dmarce1
 */

#ifndef BGRID_HPP_
#define BGRID_HPP_

#include "fwd.hpp"

namespace xtree{
template<typename T, typename Dims, int Bw = 1>
class bgrid: public grid_base<T, Dims::dim()> {
public:
	static constexpr int Ndim = Dims::dim();
	using grid_type = grid<T, Dims, Bw>;
	using index_type = typename grid_base<T, Ndim>::index_type;
private:
	const grid_type* local_ptr;
	const typename grid_type::dir_type dir;
	int size;
	std::vector<T> data;
private:
	int vector_to_index(const index_type& i) {
		if (local_ptr) {
			return grid<T, Dims, Bw>::vector_to_index(i);
		} else {
			int index = 0;
			for (int d = Ndim - 1; d >= 0; d--) {
				switch (dir[d]) {
				case -1:
					index *= Bw;
					index += i[d] + Bw - Dims::get(d);
					break;
				case 0:
					index *= Dims::get(d);
					index += i[d];
					break;
				case +1:
					index *= Bw;
					index += i[d];
					break;
				}
			}
			return index;
		}
	}
	const T& get(const index_type& i) const {
		return data[vector_to_index(i)];
	}
public:
	bgrid(const typename grid_type::dir_type& _dir, const grid_type* lptr) :
			dir(_dir), local_ptr(lptr) {
	}
	virtual ~bgrid() {
	}
	template<typename Arc>
	void load(Arc& ar, const int v) {
		local_ptr = nullptr;
		size = 1;
		for (int i = 0; i < Ndim; i++) {
			if (dir[i] == 0) {
				size *= Dims::get(i);
			} else {
				size *= Bw;
			}
		}
		data.resize(size);
		for( int i = 0; i < size; i++) {
			ar >> data[i];
		}
	}
	template<typename Arc>
	void save(Arc& ar, const int v) const {
		vector<int, Ndim> min, max;
		for (int d = 0; d < Ndim; d++) {
			switch (dir[d]) {
			case -1:
				min[d] = Dims::get(d) - Bw;
				max[d] = Dims::get(d) - 1;
				break;
			case 0:
				min[d] = 0;
				max[d] = Dims::get(d) - 1;
				break;
			case +1:
				min[d] = 0;
				max[d] = Bw - 1;
				break;
			}
		}
		for (grid_index<Ndim> i(min, max); !i.end(); ++i) {
			ar << local_ptr->get(i);
		}
	}
};

}

#endif /* BGRID_HPP_ */
