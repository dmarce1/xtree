/*
 * grid_index.hpp
 *
 *  Created on: May 29, 2014
 *      Author: dmarce1
 */

#ifndef GRID_INDEX_HPP_
#define GRID_INDEX_HPP_

#include "fwd.hpp"
namespace xtree {

template<int Ndim>
class grid_index: public std::array<int, Ndim> {
private:
	std::array<int, Ndim> min;
	std::array<int, Ndim> max;
public:
	std::array<double, Ndim> to_double() const {
		std::array<double, Ndim> d;
		for (int di = 0; di < Ndim; di++) {
			d[di] = double((*this)[di]);
		}
		return d;
	}
	grid_index(const std::array<int, Ndim>& min_, const std::array<int, Ndim>& max_) :
			std::array<int, Ndim>(min_), min(min_), max(max_) {
	}
	grid_index(const std::array<int, Ndim>& max_) :
			max(max_) {
		std::fill(std::array<int, Ndim>::begin(), std::array<int, Ndim>::end(), 0);
		std::fill(min.begin(), min.end(), 0);
	}
	grid_index(const grid_index<Ndim>& gi) :
			std::array<int, Ndim>(gi), min(gi.min), max(gi.max) {
	}
	virtual ~grid_index() {
	}
	void operator++() {
		int i = 0;
		while (((*this)[i] == max[i]) && (i != Ndim - 1)) {
			(*this)[i] = min[i];
			i++;
		}
		(*this)[i]++;
	}
	void operator++(int) {
		operator++();
	}
	bool end() {
		for (int i = 0; i < Ndim - 1; i++) {
			if ((*this)[i] != 0) {
				return false;
			}
		}
		return bool((*this)[Ndim - 1] == max[Ndim - 1] + 1);
	}
};
}

#endif /* GRID_INDEX_HPP_ */
