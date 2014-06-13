/*
 * container_math.hpp
 *
 *  Created on: Jun 4, 2014
 *      Author: dmarce1
 */

#ifndef CONTAINER_MATH_HPP_
#define CONTAINER_MATH_HPP_


template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator*(const Container<T, Size...>& a, const T& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), c.begin(), [&b](const T& this_a) {
				return this_a * b;
			});
	return c;
}

template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator+(const Container<T, Size...>& a, const T& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), c.begin(), [&b](const T& this_a) {
				return this_a + b;
			});
	return c;
}

template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator/(const Container<T, Size...>& a, const T& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), c.begin(), [&b](const T& this_a) {
				return this_a / b;
			});
	return c;
}


template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator%(const Container<T, Size...>& a, const T& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), c.begin(), [&b](const T& this_a) {
				return this_a % b;
			});
	return c;
}

template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator-(const Container<T, Size...>& a, const T& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), c.begin(), [&b](const T& this_a) {
				return this_a - b;
			});
	return c;
}

template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator*(const Container<T, Size...>& a,const Container<T, Size...>& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), [](const T& this_a, const T& this_b) {
				return this_a * this_b;
			});
	return c;
}

template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator+(const Container<T, Size...>& a,const Container<T, Size...>& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), [](const T& this_a, const T& this_b) {
				return this_a + this_b;
			});
	return c;
}

template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator/(const Container<T, Size...>& a,const Container<T, Size...>& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), [](const T& this_a, const T& this_b) {
				return this_a / this_b;
			});
	return c;
}

template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator%(const Container<T, Size...>& a,const Container<T, Size...>& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), [](const T& this_a, const T& this_b) {
				return this_a % this_b;
			});
	return c;
}

template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container>
Container<T, Size...> operator-(const Container<T, Size...>& a,const Container<T, Size...>& b) {
	Container<T, Size...> c;
	std::transform(a.begin(), a.end(), b.begin(), c.begin(), [](const T& this_a, const T& this_b) {
				return this_a - this_b;
			});
	return c;
}


#endif /* CONTAINER_MATH_HPP_ */
