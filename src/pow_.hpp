/*
 * pow_.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef POW__HPP_
#define POW__HPP_


namespace xtree {

template<int Base, int Exponent>
struct pow_ {
	static constexpr int get() {
		return Base * pow_<Base, Exponent - 1>::get();
	}
};

template<int Base>
struct pow_<Base, 0> {
	static constexpr int get() {
		return 1;
	}
};

template<int Exponent>
struct pow_<0, Exponent> {
	static constexpr int get() {
		return 1;
	}
};

}

#endif /* POW__HPP_ */
