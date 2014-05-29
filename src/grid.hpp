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
	friend class xgrid<T, Dims> ;
	friend class std::allocator<grid>;
public:
	static constexpr int Ndim = Dims::dim();
	static constexpr int Nchild = 1 << Ndim;
	static constexpr int Size = Dims::size();
	static constexpr int bw = Bw;
	using base_type = grid_base<T, Ndim>;
	using index_type = typename grid_base<T, Ndim>::index_type;
	using type = T;
	using dims_type = Dims;
private:
	std::array<T, Size> data;
	std::array<std::shared_ptr<const base_type>, pow_<3, Ndim>::value> grid_selector;
private:
	T& get(const index_type& i) {
		return data[base_type::template vector_to_index<Dims>(i)];
	}
	const T& get(const index_type& i) const {
		return data[base_type::template vector_to_index<Dims>(i)];
	}
private:
	grid() {
		for (dir_type<Ndim> i; !i.end(); i++) {
			grid_selector[i] = nullptr;
		}
		dir_type<Ndim> i;
		for (int di = 0; di < Ndim; di++) {
			i[di] = 0;
		}
	}
public:

	static std::shared_ptr<grid> create() {
		std::allocator<grid> alloc;
		std::shared_ptr < grid > ptr(new grid<T, Dims, Bw>);
		dir_type<Ndim> i;
		for (int di = 0; di < Ndim; di++) {
			i[di] = 0;
		}
		ptr->grid_selector[i] = ptr;
		return ptr;
	}
	T& operator[](const index_type& i) {
		return get(i);
	}
	const T& operator[](const index_type& global_index) const {
		dir_type<Ndim> grid_index;
		index_type local_index;
		for (int i = 0; i < Ndim; i++) {
			const div_t tmp1 = div(global_index[i] + Dims::get(i), Dims::get(i));
			grid_index[i] = tmp1.quot - 1;
			local_index[i] = tmp1.rem;
		}
		return grid_selector[grid_index]->get(local_index);
	}
	void set_neighbor(const std::shared_ptr<const base_type>& n, dir_type<Ndim> dir) {
		grid_selector[dir] = n;
	}
	virtual ~grid() {
	}
	void prolong(const xgrid<T, Dims>& x) {
		for (array_index<Ndim> i(Dims::to_vector() - 1); !i.end(); i++) {
			(*this)[i] = x[i / 2];
		}
	}
	vector<T, Size / Nchild> restrict_to_stream() const {
		constexpr T factor = T(1) / T(Nchild);
		using this_dims = int_seq_over2<Dims>;
		vector<T, Size / Nchild> stream(T(0));
		const vector<int, Ndim> this_max = this_dims::to_vector() - 1;
		for (array_index<Ndim> ci(this_max); !ci.end(); ci++) {
			for (child_index_type<Ndim> xi; !xi.end; xi++) {
				const vector<int, Ndim> fi = ci * 2 + xi.to_vector();
				stream[base_type::template vector_to_index<this_dims>(ci)] += (*this)[fi] * factor;
			}
		}
		return stream;
	}
	void restrict_from_stream(const child_index_type<Ndim>& xi, const vector<T, Size / Nchild>& stream) {
		const vector<int, Ndim> this_min = (xi.to_vector() * Dims::to_vector()) / 2;
		const vector<int, Ndim> this_max = this_min + (Dims::to_vector() / 2) - 1;
		int j = 0;
		for (array_index<Ndim> ci(this_min, this_max); !ci.end(); ci++) {
			(*this)[ci] = stream[j];
			j++;
		}
	}
};

}

#endif /* ARRAY_HPP_ */
