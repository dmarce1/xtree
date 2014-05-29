/*
 * agrid.hpp
 *
 *  Created on: May 29, 2014
 *      Author: dmarce1
 */

#ifndef AGRID_HPP_
#define AGRID_HPP_

#include "fwd.hpp"

namespace xtree {
template<typename T, typename Dims, int Bw = 1>
class agrid: public grid_base<T, Dims::dim()> {
public:
	static constexpr int Ndim = Dims::dim();
	static constexpr int Abw = (Bw - 1) / 2 + 1;
	using grid_type = grid<T, Dims, Abw>;
	using this_dims = int_seq_over2<Dims>;
	using index_type = typename grid_base<T, Ndim>::index_type;
private:
	const grid_type* local_ptr;
	const dir_type<Ndim> dir;
	const child_index_type<Ndim> xi;
	vector<int, Ndim> offset;
	int size;
	std::vector<T> data;
private:
	const T& get(const index_type& i) const {
		if (local_ptr) {
			return (*local_ptr)[offset + i / 2];
		} else {
			int index = 0;
			for (int d = Ndim - 1; d >= 0; d--) {
				switch (dir[d]) {
				case 0:
					index *= this_dims::get(d);
					index += i[d] / 2;
					break;
				case -1:
					index *= Abw;
					index += i[d] / 2 + Abw - this_dims::get(d);
					break;
				case +1:
					index *= Abw;
					index += i[d] / 2;
					break;
				default:
					assert(false);
					break;
				}
			}
			return data[index];
		}
	}
public:
	agrid() {
	}
	agrid(const child_index_type<Ndim>& _xi, const dir_type<Ndim>& _dir, const grid_type* lptr) :
			xi(_xi), dir(_dir), local_ptr(lptr) {
		for (int i = 0; i < Ndim; i++) {
			if (dir[i] == 0) {
				offset[i] = xi[i] * this_dims::get(i);
			} else {
				offset[i] = 0;
			}
		}
	}
	virtual ~agrid() {
	}
	template<typename Arc>
	void load(Arc& ar, const int v) {
		local_ptr = nullptr;
		size = 1;
		for (int i = 0; i < Ndim; i++) {
			if (dir[i] == 0) {
				size *= this_dims::get(i);
			} else {
				size *= Abw;
			}
		}
		data.resize(size);
		data.load(ar, v);
	}
	template<typename Arc>
	void save(Arc& ar, const int v) const {
		if (local_ptr) {
			vector<int, Ndim> min, max;
			for (int d = 0; d < Ndim; d++) {
				switch (dir[d]) {
				case 0:
					min[d] = offset[d];
					max[d] = this_dims::get(d) + offset[d] - 1;
					break;
				case -1:
					min[d] = this_dims::get(d) - Abw;
					max[d] = this_dims::get(d) - 1;
					break;
				case +1:
					min[d] = 0;
					max[d] = Abw - 1;
					break;
				default:
					assert(false);
					break;
				}
			}
			for (grid_index<Ndim> i(min, max); !i.end(); ++i) {
				ar << local_ptr->get(i);
			}
		} else {
			data.save(ar, v);
		}
	}
};

}

#endif /* AGRID_HPP_ */