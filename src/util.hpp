/*
 * int_seq.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef INT_SEQ_HPP_
#define INT_SEQ_HPP_

namespace xtree {

template<typename U, typename T, int Ndim>
using get_type = T (U::*)(location<Ndim>);

template<typename U, typename T, int Ndim>
using set_type = void (U::*)(location<Ndim>, T);

enum op_type {
	REBRANCH, ASCEND, DESCEND, EXCHANGE
};


template<int N, int ...Params>
struct int_seq {
	static constexpr int dim() {
		return 1 + sizeof...(Params);
	}
	static constexpr int get(const int i) {
		return i == 0 ? N : int_seq<Params...>::get(i - 1);
	}
	static constexpr int size() {
		return N * int_seq<Params...>::size();
	}
};

template<int N>
struct int_seq<N> {
	static constexpr int dim() {
		return 1;
	}
	static constexpr int get(const int i) {
		return N;
	}
	static constexpr int size() {
		return N;
	}
};



template<int Base, int Exponent>
struct pow_ {
	static constexpr int value = Base * pow_<Base, Exponent - 1>::value;
};

template<int Base>
struct pow_<Base, 0> {
	static constexpr int value = 1;
};

}

#endif /* INT_SEQ_HPP_ */
