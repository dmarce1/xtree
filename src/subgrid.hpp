/*
 * subgrid.hpp
 *
 *  Created on: May 29, 2014
 *      Author: dmarce1
 */

#ifndef SUBGRID3_HPP_
#define SUBGRID3_HPP_

#include "fwd.hpp"

namespace xtree {

template<typename T, typename Dims>
class subgrid: public grid_base<T, Dims::dim()> {
public:
	static constexpr int Ndim = Dims::dim();
	using base_type = grid_base<T, Ndim>;
	using index_type = typename grid_base<T, Ndim>::index_type;
private:
	static constexpr int Size = Dims::size();
	std::array<T, Size> data;
private:
	T& get(const index_type& i) {
		return data[base_type::vector_to_index(i)];
	}
	const T& get(const index_type& i) const {
		return data[base_type::vector_to_index(i)];
	}
public:
	subgrid() {
	}
	T& operator[](const index_type& i) {
		return get(i);
	}
	const T& operator[](const index_type& global_index) const {
		return get(i);
	}
	virtual ~subgrid() {
	}
	template< typename Arc >
	void serialize( Arc& ar, const int v ) {
		data.serialize(ar,v);
	}
};

}

#endif /* SUBGRID_HPP_ */
