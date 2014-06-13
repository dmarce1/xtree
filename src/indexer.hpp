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
class indexer: public std::array<int, Ndim> {
private:
	bool is_end;
private:
	int abs_value(int i) const {
		return (*this)[i] - Origin;
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
			(*this)[j] = d.rem;
			j = d.quot;
		}
		return *this;
	}
	template<typename Archive>
	void serialize(Archive& ar, const int v) {
		ar & is_end;
		for (int i = 0; i < Ndim; i++) {
			ar & (*this)[i];
		}
	}
	int operator[](int i) const {
		return (*this)[i];
	}
	int& operator[](int i) {
		return (*this)[i];
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
				is_end = true;
				return;
			}
		}
		(*this)[i]++;
	}
	void begin() {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] = Origin;
		}
		is_end = false;
	}
	bool end() {
		return is_end;
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
	void set_zero() {
		for (int i = 0; i < Ndim; i++) {
			std::array<int, Ndim>::operator[](i) = 0;
		}
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
