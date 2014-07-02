/*
 * indexer.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef INDEXER_HPP_
#define INDEXER_HPP_

namespace xtree {

template<int Ndim, int Size, int Origin = 0>
class indexer {
protected:
	std::array<int, Ndim> values;
	bool __is_end;
	int abs_value(int i) const {
		return values[i] - Origin;
	}
public:
	indexer(const std::array<int, Ndim>& is) {
		*this = is;
	}
	indexer() {
		begin();
	}
	indexer& operator=(int i) {
		for (int j = 0; j < Ndim; j++) {
			auto d = div(i, Size);
			(*this)[j] = d.rem + Origin;
			i = d.quot;
		}
		return *this;
	}
	template<typename Archive>
	void serialize(Archive& ar, const int v) {
		ar & __is_end;
		for (int i = 0; i < Ndim; i++) {
			ar & (*this)[i];
		}
	}
	int operator[](int i) const {
		return values[i];
	}
	int& operator[](int i) {
		return values[i];
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
			(*this)[i] = Origin;
			i++;
			if (i == Ndim) {
				__is_end = true;
				return;
			}
		}
		(*this)[i]++;
	}
	void begin() {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] = Origin;
		}
		__is_end = false;
	}
	bool is_end() {
		return __is_end;
	}
	indexer& flip(int i) {
		(*this)[i] = -(*this)[i] + 2 * Origin + Size - 1;
		return *this;
	}
	indexer& flip() {
		for (int i = 0; i < Ndim; i++) {
			flip(i);
		}
		return *this;
	}
	std::array<int, Ndim> to_vector() const {
		std::array<int, Ndim> v;
		for (int di = 0; di < Ndim; di++) {
			v[di] = (*this)[di];
		}
		return v;
	}
	indexer& set_zero() {
		for (int i = 0; i < Ndim; i++) {
			values[i] = 0;
		}
		return *this;
	}
	bool is_zero() {
		bool rc = true;
		for (int i = 0; i < Ndim; i++) {
			if (std::array<int, Ndim>::operator[](i) != 0) {
				rc = false;
				break;
			}
		}
		return rc;
	}
};

}
#endif /* INDEXER_HPP_ */
