/*
 * valarray.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: dmarce1
 */

#ifndef VALAR24RAY_HPP_
#define VALAR24RAY_HPP_

#include <valarray>


std::gslice get_slice(const std::valarray<std::size_t>& dims, const std::valarray<std::size_t>& mins, const std::valarray<std::size_t>& maxes) {
	std::valarray<std::size_t> strides;
	const std::size_t ndim = dims.size();
	std::size_t start;
	const std::valarray<std::size_t> these_dims(maxes - mins);
	strides[0] = 1;
	start = mins[0];
	for (std::size_t i = 1; i < ndim; i++) {
		strides[i] = strides[i - 1] * dims[i - 1];
		start += strides[i] * mins[i];
	}
	return std::gslice(start, these_dims, strides);
}

std::gslice get_xtant_slice(const std::valarray<std::size_t>& dims, std::size_t ci) {
	std::valarray<std::size_t> mins(dims.size());
	std::valarray<std::size_t> maxes(dims.size());
	for (std::size_t d = 0; d < dims.size(); d++) {
		if ((ci >> d) & 1) {
			mins[d] = dims[d] / std::size_t(2);
			maxes[d] = dims[d];
		} else {
			mins[d] = std::size_t(0);
			maxes[d] = dims[d] / std::size_t(2);
		}
	}
	return get_slice(dims, mins, maxes);

}

std::gslice get_restrict_slice(const std::valarray<std::size_t>& dims, std::size_t ci) {
	std::valarray<std::size_t> strides(dims.size());
	std::valarray<std::size_t> these_dims = dims / std::size_t(2);
	strides[0] = 1;
	for (std::size_t d = 1; d < dims.size(); d++) {
		strides[d] = dims[d - 1] * strides[d - 1];
	}
	std::size_t start = 0;
	for (std::size_t d = 0; d < dims.size(); d++) {
		start += strides[d] * ((ci >> d) & 1);
	}
	strides *= std::size_t(2);
	return std::gslice(start, these_dims, strides);

}

template<typename T>
std::valarray<T> get_prolong_array(const std::valarray<T>& a, const std::valarray<std::size_t>& dims) {
	const auto nchild = 1 << dims.size();
	std::valarray<T> b(nchild);
	for (std::size_t ci = 0; ci < nchild; ci++) {
		b[get_restrict_slice(dims, ci)] = a;
	}
	return b;
}

template<typename T>
T product(std::valarray<T> dims) {
	std::size_t p = dims[0];
	for (std::size_t i = 1; i < dims.size(); i++) {
		p *= dims[i];
	}
	return p;
}

std::valarray<std::valarray<double>> create_position_array(std::valarray<std::size_t> dims) {
	std::valarray<std::valarray<double>> X;
	std::valarray<std::size_t> strides(dims.size());
	X.resize(product(dims));
	std::valarray<std::valarray<double>> unit(dims.size());
	for (std::size_t i = 0; i < dims.size(); i++) {
		unit[i].resize(dims.size(), 0.0);
		unit[i][i] = 1.0;
	}
	for (std::size_t i = 0; i < product(dims); i++) {
		X[i].resize(dims.size(), 0.0);
	}
	strides[0] = 1;
	for (std::size_t di = 1; di < dims.size(); di++) {
		strides[di] = strides[di - 1] * dims[di - 1];
	}
	for (std::size_t di = 0; di < dims.size(); di++) {
		std::valarray<std::size_t> these_dims(dims);
		dims[di] = 1;
		for (std::size_t i = 0; i < dims[di]; i++) {
			const auto start = i * strides[di];
			const auto slice = std::gslice(start, dims, strides);
			std::valarray<std::valarray<double>> x(product(dims) / dims[di]);
			x = std::valarray<double>(unit[di] * (double(i) + 0.5));
			X[slice] += x;
		}
	}
	return X;
}



#endif /* VALARRAY_HPP_ */
