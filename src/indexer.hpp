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
private:
	std::array<int, Ndim> value;
	bool is_end;
private:
	int abs_value(int i) const {
		return value[i] - Origin;
	}
public:
	indexer() {
		begin();
	}
	indexer& operator=(int i) {
		for (int j = 0; j < Ndim; j++) {
			auto d = div(i, Size);
			value[j] = d.rem;
			j = d.quot;
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
	int operator[](int i) const {
		return value[i];
	}
	int& operator[](int i) {
		return value[i];
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
	void begin() {
		for (int i = 0; i < Ndim; i++) {
			value[i] = Origin;
		}
		is_end = false;
	}
	bool end() {
		return is_end;
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
	std::array<int, Ndim> to_vector() const {
		std::array<int, Ndim> v;
		for (int di = 0; di < Ndim; di++) {
			v[di] = value[di];
		}
		return v;
	}
};

}
#endif /* INDEXER_HPP_ */
