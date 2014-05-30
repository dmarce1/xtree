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
using get_type = T (U::*)(const location<Ndim>&);

template<typename U, typename T, int Ndim>
using set_type = void (U::*)(const location<Ndim>&, T);

template<int Ndim>
using child_index_type = indexer<Ndim, 2, 0>;

template<int Ndim>
using dir_type = indexer<Ndim, 3, -1>;

enum bound_type {
	DECOMP, AMR, PHYS
};

enum op_type {
	REBRANCH, ASCEND, DESCEND, EXCHANGE, AMR_ASCEND
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
	static vector<int, dim()> to_vector() {
		vector<int, dim()> v;
		for (int i = 0; i < dim(); i++) {
			v[i] = get(i);
		}
		return v;
	}
};

template<typename Sequence>
struct int_seq_over2 {
	static constexpr int dim() {
		return Sequence::dim();
	}
	static constexpr int get(const int i) {
		return Sequence::get(i) / 2;
	}
	static constexpr int size() {
		return Sequence::size() / (1 << dim());
	}
	static vector<int, dim()> to_vector() {
		vector<int, dim()> v;
		for (int i = 0; i < dim(); i++) {
			v[i] = get(i);
		}
		return v;
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
	static vector<int, dim()> to_vector() {
		return vector<int, dim()>(N);
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
