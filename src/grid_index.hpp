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
class grid_index: public vector<int, Ndim> {
private:
	const vector<int, Ndim> min;
	const vector<int, Ndim> max;
public:
	grid_index(const vector<int, Ndim>& min_, const vector<int, Ndim>& max_) :
			vector<int, Ndim>(min_), min(min_), max(max_) {
	}
	grid_index(const vector<int, Ndim>& max_) :
			vector<int, Ndim>(0), min(0), max(max_) {
	}
	void operator++() {
		int i = 0;
		while ((*this)[i] == max[i]) {
			(*this)[i] = min[i];
			i++;
			assert(i != Ndim);
		}
		(*this)[i]++;
	}

	void operator++(int) {
		operator++();
	}
	bool end() {
		for (int i = 0; i < Ndim; i++) {
			if ((*this)[i] != max[i] + i / (Ndim - 1)) {
				return false;
			}
		}
		return true;
	}
	virtual ~grid_index() {
	}
};
}

#endif /* GRID_INDEX_HPP_ */
