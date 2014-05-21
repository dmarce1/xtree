/*
 * indexer.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef INDEXER_HPP_
#define INDEXER_HPP_

#include "xtree.hpp"
#include "int_seq.hpp"

namespace xtree {
template<typename Dims, class origin = int_seq_const<0, Dims::dim()> >
class indexer {
private:
	static constexpr int Ndim = Dims::dim();
	int value[Ndim];
	bool is_end;
	int abs_value(int i) const {
		return value[i] - origin::get(i);
	}
public:
	int operator[](int i) const {
		return value[i];
	}
	int& operator[](int i) {
		return value[i];
	}
	indexer() {
		begin();
	}
	void begin() {
		for (int i = 0; i < Ndim; i++) {
			value[i] = origin::get(i);
		}
		is_end = false;
	}
	bool end() {
		return is_end;
	}
	operator int() const {
		int j = abs_value(Ndim - 1);
		for (int i = Ndim - 2; i >= 0; i--) {
			j *= Dims::get(i);
			j += abs_value(i);
		}
		return j;
	}
	void operator++(int) {
		int i = 0;
		while (abs_value(i) == Dims::get(i) - 1) {
			value[i] = origin::get(i);
			i++;
			if (i == Dims::dim()) {
				is_end = true;
				return;
			}
		}
		value[i]++;
	}
	indexer& flip(int i) {
		value[i] = -value[i] + 2 * origin::get(i) + Dims::get(i) - 1;
		return *this;
	}
	indexer& flip() {
		for (int i = 0; i < Ndim; i++) {
			flip(i);
		}
		return *this;
	}
	template<typename Archive>
	void serialize(Archive& ar, const int v) {
		ar & is_end;
		for (int i = 0; i < Ndim; i++) {
			ar & value[i];
		}
	}
};
}
#endif /* INDEXER_HPP_ */
