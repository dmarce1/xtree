/*
 * int_seq.hpp
 *
 *  Created on: May 19, 2014
 *      Author: dmarce1
 */

#ifndef INT_SEQ_HPP_
#define INT_SEQ_HPP_

#include <thread>

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
	static constexpr int value = N;
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
	static std::array<int, dim()> to_vector() {
		std::array<int, dim()> v;
		for (int i = 0; i < dim(); i++) {
			v[i] = get(i);
		}
		return v;
	}
};

template<typename Sequence, int Incdim>
struct int_seq_plus_one {
	static constexpr int dim() {
		return Sequence::dim();
	}
	static constexpr int get(const int i) {
		return Sequence::get(i) + ((i == Incdim) ? 1 : 0);
	}
	static constexpr int size() {
		return Sequence::size() + Sequence::size() / Sequence::get(Incdim);
	}
	static std::array<int, dim()> to_vector() {
		std::array<int, dim()> v;
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
	static std::array<int, dim()> to_vector() {
		std::array<int, dim()> v;
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
	static std::array<int, dim()> to_vector() {
		return std::array<int, dim()>(N);
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
	static constexpr std::size_t value = factorial_<N - 1>::value * N;
};

template<>
struct factorial_<0> {
	static constexpr std::size_t value = 1;
};
}

template<class R, class ...Args>
R exec_on_separate_thread(R (*fptr)(Args...), Args ...args) {
	R data;
	std::thread([&](Args*...args) {
		data = (*fptr)(*args...);
	}, &args...).join();
	return data;

}

template<typename T>
class future_wrapper {
	mutable std::shared_ptr<hpx::future<T>> ptr;
public:
	future_wrapper() :
			ptr(std::make_shared<hpx::future<T>>()) {
	}
	~future_wrapper() = default;
	future_wrapper(const future_wrapper&) = default;
	future_wrapper& operator=(const future_wrapper&) = default;
	future_wrapper(hpx::future<T>&& f) :
			ptr(std::make_shared < hpx::future < T >> (std::move(f))) {
	}
	future_wrapper& operator=(hpx::future<T>&& f) {
		*ptr = std::move(f);
		return *this;
	}
	T get() {
		return ptr->get();
	}
	hpx::shared_future<T> share() {
		return ptr->share();
	}

	template<typename Arc>
	void load(Arc& ar, const unsigned v) {
		T data;
		ar & data;
		*ptr = hpx::make_ready_future < T > (std::move(data));
	}

	template<typename Arc>
	void save(Arc& ar, const unsigned v) const {
		T data = ptr->get();
		ar & data;
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER();
};

#endif /* INT_SEQ_HPP_ */
