/*
 * indexer.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef INDEXER_HPP_
#define INDEXER_HPP_

namespace xtree {

template<int Ndim, int Size, int Origin=0 >
class indexer {
private:
	int value[Ndim];
	bool is_end;
	int abs_value(int i) const {
		return value[i] - Origin;
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
			value[i] = Origin;
		}
		is_end = false;
	}
	bool end() {
		return is_end;
	}
	operator int() const {
		int j = abs_value(Ndim - 1);
		for (int i = Ndim - 2; i >= 0; i--) {
			j *= Size;
			j += abs_value(i);
		}
		return j;
	}
	void operator++(int) {
		int i = 0;
		while (abs_value(i) == Size - 1) {
			value[i] = Origin;
			i++;
			if (i == Ndim) {
				is_end = true;
				return;
			}
		}
		value[i]++;
	}
	indexer& flip(int i) {
		value[i] = -value[i] + 2 * Origin + Size - 1;
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
