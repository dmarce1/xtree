/*
 * int_seq.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef INT_SEQ_HPP_
#define INT_SEQ_HPP_

namespace xtree {

struct nullclass {
};

//template<typename U, typename T, int Ndim>
//using get_type = T (U::*)(const location<Ndim>&);

//template<typename U, typename T, int Ndim>
//using set_type = void (U::*)(const location<Ndim>&, T);
enum bound_type {
  DECOMP, AMR, PHYS
};

enum op_type {
  REBRANCH, ASCEND, DESCEND, EXCHANGE, AMR_ASCEND, LOCAL
};

template<int N>
struct int2type {
  static const int value = N;
};

template<int N, int ...Params>
struct int_seq {
    static const int dim = 1 + sizeof...(Params);
    
//    static const int dim() {
//		return 1 + sizeof...(Params);
//	}
  static const int get(const int i) {
    return i == 0 ? N : int_seq<Params...>::get(i - 1);
  }
  static const int size() {
    return N * int_seq<Params...>::size();
  }
  static std::array<int, dim> to_vector() {
    std::array<int, dim> v;
    for (int i = 0; i < dim; i++) {
      v[i] = get(i);
    }
    return v;
  }
};

template<typename Sequence, int Incdim>
struct int_seq_plus_one {
    static const int dim = Sequence::dim;
    
//    static const int dim() {
//		return Sequence::dim();
//	}
  static const int get(const int i) {
    return Sequence::get(i) + ((i == Incdim) ? 1 : 0);
  }
  static const int size() {
    return Sequence::size() + Sequence::size() / Sequence::get(Incdim);
  }
  static std::array<int, dim> to_vector() {
    std::array<int, dim> v;
    for (int i = 0; i < dim; i++) {
      v[i] = get(i);
    }
    return v;
  }
};

template<typename Sequence>
struct int_seq_over2 {
    static const int dim = Sequence::dim;

//    static const int dim() {
//		return Sequence::dim();
//	}
  static const int get(const int i) {
    return Sequence::get(i) / 2;
  }
  static const int size() {
    return Sequence::size() / (1 << dim());
  }
  static std::array<int, dim> to_vector() {
    std::array<int, dim> v;
    for (int i = 0; i < dim; i++) {
      v[i] = get(i);
    }
    return v;
  }
};

template<int N>
struct int_seq<N> {
    static const int dim = 1;

//    static const int dim() {
//		return 1;
//	}
  static const int get(const int i) {
    return N;
  }
  static const int size() {
    return N;
  }
  static std::array<int, dim> to_vector() {
    return std::array<int, dim>(N);
  }
};

template<int Base, int Exponent>
struct pow_ {
  static const int value = Base * pow_<Base, Exponent - 1>::value;
};

template<int Base>
struct pow_<Base, 0> {
  static const int value = 1;
};

template<typename T>
T factorial(T i) {
  if (i != 0) {
    return i * factorial(i - 1);
  } else {
    return 1;
  }
}

template<std::size_t N>
struct factorial_ {
  static const std::size_t value = factorial_<N - 1>::value * N;
};

template<>
struct factorial_<0> {
  static const std::size_t value = 1;
};

}

#endif /* INT_SEQ_HPP_ */
