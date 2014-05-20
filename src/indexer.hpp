/*
 * indexer.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef INDEXER_HPP_
#define INDEXER_HPP_

#include "xtree.hpp"

namespace xtree {
template<typename Dims>
class indexer {
private:
	static constexpr int Ndim = Dims::dim();
	int value[Ndim];
public:
	indexer() {
		begin();
	}
	void begin() {
		for (int i = 0; i < Ndim; i++) {
			value[i] = 0;
		}
	}
	bool end() {
		for (int i = 0; i < Ndim; i++) {
			if (value[i] < Dims::get(i)) {
				return false;
			}
		}
		return true;
	}
	operator int() const {
		int j = value[Ndim - 1];
		for (int i = Ndim - 2; i >= 0; i--) {
			j *= Dims::get(i);
			j += value[i];
		}
		return j;
	}
	void operator++(int) {
		int i = 0;
		while (value[i] == Dims::get(i) - 1) {
			value[i] = 0;
			i++;
		}
		value[i]++;
	}
	indexer& flip(int i) {
		value[i] = Dims::get(i) - 1 - value[i];
		return *this;
	}
	template<typename Archive>
	void serialize(Archive& ar, const int v) {
		for (int i = 0; i < Ndim; i++) {
			ar & value[i];
		}
	}
};
}
#endif /* INDEXER_HPP_ */
