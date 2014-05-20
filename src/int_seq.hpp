/*
 * int_seq.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef INT_SEQ_HPP_
#define INT_SEQ_HPP_

#include "xtree.hpp"

namespace xtree {

template<int N, int ...Params>
struct int_seq {
	static constexpr int dim() {
		return 1 + sizeof...(Params);
	}
	static constexpr int get(int i) {
		return i == 0 ? N : int_seq<Params...>::get(i - 1);
	}
};

template<int N>
struct int_seq<N> {
	static constexpr int dim() {
		return 1;
	}
	static constexpr int get(int i) {
		return N;
	}
};

template<int N, int Ndim>
struct int_seq_const {
	static constexpr int dim() {
		return Ndim;
	}
	static constexpr int get(int i) {
		return N;
	}
};


}

#endif /* INT_SEQ_HPP_ */
