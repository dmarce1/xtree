/*
 * expansion.hpp
 *
 *  Created on: Jun 19, 2014
 *      Author: dmarce1
 */

#ifndef EXPANSION_HPP_
#define EXPANSION_HPP_

#include <valarray>

namespace xtree {

template<std::size_t P>
struct expansion: public std::valarray<complex> {
private:
	std::valarray<complex>& base_ref;
public:
	operator std::valarray<complex>() {
		return base_ref;
	}
	operator std::valarray<complex>&() {
		return base_ref;
	}
	expansion() :
			std::valarray<complex>(P * (P + 1) / 2), base_ref(*static_cast<std::valarray<complex>*>(this)) {
	}
	template<typename T>
	void serialize(T& arc, const unsigned v) {
		boost::serialization::serialize(arc, *static_cast<std::valarray<complex>*>(this), v);
	}
	expansion(const expansion& other) :
			std::valarray<complex>(P * (P + 1) / 2), base_ref(*static_cast<std::valarray<complex>*>(this)) {
		base_ref = other;
	}
	expansion(expansion&& other) :
			std::valarray<complex>(P * (P + 1) / 2), base_ref(*static_cast<std::valarray<complex>*>(this)) {
		base_ref = std::move(other);
	}
	expansion& operator=(const expansion& other) {
		base_ref = other;
		return *this;
	}
	expansion& operator=(expansion&& other) {
		base_ref = std::move(other);
		return *this;
	}
	expansion& operator=(real zero) {
		std::fill(std::begin(*this), std::end(*this), zero);
		return *this;
	}
	expansion(const std::valarray<complex>& other) :
			std::valarray<complex>(P * (P + 1) / 2), base_ref(*static_cast<std::valarray<complex>*>(this)) {
		std::copy(std::begin(other), std::end(other), std::begin(*this));
	}
};

}

#endif /* EXPANSION_HPP_ */
