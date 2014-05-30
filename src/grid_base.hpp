/*
 * grid_base.hpp
 *
 *  Created on: May 29, 2014
 *      Author: dmarce1
 */

#ifndef SUBGRID_BASE_HPP_
#define SUBGRID_BASE_HPP_

#include "fwd.hpp"
namespace xtree {

template<typename T, int Ndim>
class grid_base {
public:
	using index_type = vector<int, Ndim>;
	virtual ~grid_base() {
	}
protected:
	template< typename Dims >
	static int vector_to_index(const index_type& i) {
		int index = i[Ndim - 1];
		for (int d = Ndim - 2; d >= 0; d--) {
			index *= Dims::get(d);
			index += i[d];
		}
		return index;
	}

public:
	virtual const T& get(const index_type& i) const = 0;
};

}

#endif /* SUBGRID_BASE_HPP_ */
