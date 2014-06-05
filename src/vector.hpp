/*
 * vector.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */
/*
#ifndef VECTOR21_HPP_
#define VECTOR21_HPP

namespace xtree {

template<typename T, int ...N>
class vector: public std::array<T, N...> {
	static constexpr int Ndim = boost::mpl::int_<N...>::value;
public:
	vector() {
	}
	T product() const {
		T p = (*this)[0];
		for (int i = 1; i < Ndim; i++) {
			p *= (*this)[i];
		}
		return p;
	}
	vector(const T& a) {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] = a;
		}
	}
	vector(const std::vector<T>& v) {
		*this = v;
	}
	vector& operator<<=(const T& a) {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] <<= a;
		}
		return *this;
	}
	vector& operator>>=(const T& a) {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] >>= a;
		}
		return *this;
	}
	vector& operator+=(const vector& a) {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] += a[i];
		}
		return *this;
	}
	vector& operator-=(const vector& a) {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] -= a[i];
		}
		return *this;
	}
	vector& operator*=(const T& a) {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] *= a;
		}
		return *this;
	}
	vector& operator/=(const T& a) {
		for (int i = 0; i < Ndim; i++) {
			(*this)[i] /= a;
		}
		return *this;
	}
	vector operator+(const vector& a) const {
		vector<T, Ndim> b = *this;
		b += a;
		return b;
	}
	vector operator-(const vector& a) const {
		vector<T, Ndim> b = *this;
		b -= a;
		return b;
	}
	vector operator*(const T& a) const {
		vector<T, Ndim> b = *this;
		b *= a;
		return b;
	}
	vector operator*(const vector<T, Ndim>& a) const {
		vector<T, Ndim> b;
		for (int i = 0; i < Ndim; i++) {
			b[i] = (*this)[i] * a[i];
		}
		return b;
	}
	vector operator/(const T& a) const {
		vector<T, Ndim> b = *this;
		b /= a;
		return b;
	}
	bool operator==(const vector& b) const {
		for (int i = 0; i < Ndim; i++) {
			if ((*this)[i] != b[i]) {
				return false;
			}
		}
		return true;
	}
	bool operator!=(const vector& b) const {
		return !(*this == b);
	}
	template<typename Arc>
	void serialize(Arc& ar, const int v) {
		for (int i = 0; i < Ndim; i++) {
			ar & (*this)[i];
		}
	}
	template<typename Arc>
	void serialize(Arc& ar, const int v) const {
		for (int i = 0; i < Ndim; i++) {
			ar & (*this)[i];
		}
	}


	template<typename U>
	operator vector<U,N...>() {
		vector<U, N...>  v;
		for (int i = 0; i < Ndim; i++) {
			v[i] = U((*this)[i]);
		}
		return v;
	}

};

template<typename T>
class vector<T> : public std::vector<T> {
public:
	vector(std::vector<T> vec) :
			std::vector<T>(vec) {
	}
	vector() :
			std::vector<T>() {
	}
	vector(const int sz) :
			std::vector<T>(sz) {
	}
	vector& operator=(const T& a) {
		for (int i = 0; i < this->size(); i++) {
			(*this)[i] = a;
		}
		return *this;
	}
	vector& operator+=(const vector& a) {
		for (int i = 0; i < this->size(); i++) {
			(*this)[i] += a[i];
		}
		return *this;
	}
	vector& operator-=(const vector& a) {
		for (int i = 0; i < this->size(); i++) {
			(*this)[i] -= a[i];
		}
		return *this;
	}
	vector& operator*=(const T& a) {
		for (int i = 0; i < this->size(); i++) {
			(*this)[i] *= a;
		}
		return *this;
	}
	vector& operator/=(const T& a) {
		for (int i = 0; i < this->size(); i++) {
			(*this)[i] /= a;
		}
		return *this;
	}
	vector operator+(const vector& a) const {
		vector<T> b(*this);
		b += a;
		return b;
	}
	vector operator-(const vector& a) const {
		vector<T> b(*this);
		b -= a;
		return b;
	}
	vector operator*(const T& a) const {
		vector<T> b(*this);
		b *= a;
		return b;
	}
	vector operator/(const T& a) const {
		vector<T> b(*this);
		b /= a;
		return b;
	}
	bool operator==(const vector& b) const {
		for (int i = 0; i < this->size(); i++) {
			if ((*this)[i] != b[i]) {
				return false;
			}
		}
		return true;
	}
	bool operator!=(const vector& b) const {
		return !(*this == b);
	}
	T product() const {
		T p = (*this)[0];
		for (int i = 1; i < this->size(); i++) {
			p *= (*this)[i];
		}
		return p;
	}
};


}
#endif

*/
