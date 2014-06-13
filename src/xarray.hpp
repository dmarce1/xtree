/*
 * xarray.hpp
 *
 *  Created on: Jun 11, 2014
 *      Author: dmarce1
 */

#ifndef XARRAY_HPP_
#define XARRAY_HPP_

#include <boost/serialization/valarray.hpp>

namespace xtree {

template<std::size_t Ndim>
class xarray_boundary;

template<std::size_t Ndim>
class xarray: public std::valarray<double> {
	friend class xarray_boundary<Ndim> ;
private:
	std::vector<std::size_t> dims;
	std::vector<std::size_t> strides;
	std::size_t size;

	std::valarray<double>& get_base() {
		return *(static_cast<std::valarray<double>*>(this));
	}

public:
	xarray() = default;
	xarray(const xarray&) = default;
	xarray(xarray&&) = default;
	virtual ~xarray() = default;
	xarray& operator=(const xarray&) = default;
	xarray& operator=(xarray&&) = default;

	template<typename Container>
	void resize(const Container& _dims) {
		dims.resize(Ndim + 1);
		strides.resize(Ndim + 1);
		std::copy(_dims.begin(), _dims.end(), dims.begin());
		size = dims[0];
		strides[0] = 1;
		for (std::size_t i = 1; i < Ndim + 1; i++) {
			strides[i] = strides[i - 1] * dims[i - 1];
			size *= dims[i];
		}
		get_base().resize(size);
	}

	template<typename Arc>
	void serialize(Arc& ar, const unsigned int v) {
		ar & dims;
		ar & strides;
		ar & size;
		boost::serialization::serialize(ar, get_base(), v);
	}

	std::gslice get_slice(const dir_type<Ndim>& dir, std::size_t bw) const {
		std::vector<std::size_t> these_dims(Ndim + 1);
		these_dims[Ndim] = dims[Ndim];
		const auto& these_strides = strides;
		std::size_t start = 0;
		for (std::size_t i = 0; i < Ndim; i++) {
			if (dir[i] == 0) {
				these_dims[i + 1] = dims[i + 1] - 2 * bw;
				start += bw * these_strides[i];
			} else {
				these_dims[i + 1] = bw;
			}
			if (dir[i] > 0) {
				start += these_strides[i + 1] - 2 * bw * these_strides[i];
			}
		}
		auto slice = std::gslice(start, these_dims, these_strides);
		xarray<Ndim> b;
		b.resize(these_dims);
		b.get_base() = get_base()[slice];
		return b;
	}
}
;

template<std::size_t Ndim>
class xarray_boundary {

private:
	dir_type<Ndim> dir;
	std::shared_ptr<xarray<Ndim>> array_ptr;
	std::size_t bw;
	bool isvirtual;

public:
	xarray_boundary() {
		isvirtual = false;
	}

	~xarray_boundary() = default;

	void set(const dir_type<Ndim>& _dir, std::size_t _bw, xarray<Ndim>& array) {
		array_ptr = std::shared_ptr < xarray < Ndim >> (&array, [](xarray<Ndim>*) {
		});
		dir = _dir;
		bw = _bw;
	}

	xarray<Ndim> get() const {
		if (isvirtual) {
			return std::move(*array_ptr);
		} else {
			std::vector<std::size_t> these_dims(Ndim + 1);
			these_dims[Ndim] = array_ptr->dims[Ndim];
			const auto& these_strides = array_ptr->strides;
			std::size_t start = 0;
			for (std::size_t i = 0; i < Ndim; i++) {
				if (dir[i] != 0) {
					these_dims[i + 1] = bw;
				} else {
					these_dims[i + 1] = array_ptr->dims[i + 1];
				}
				if (dir[i] > 0) {
					start += these_strides[i + 1] - 2 * bw * these_strides[i];
				}
			}
			auto slice = std::gslice(start, these_dims, these_strides);
			xarray<Ndim> b;
			b.resize(these_dims);
			b.get_base() = (*array_ptr)[slice];
			return b;
		}
	}

	template<typename Arc>
	void save(Arc& ar, const unsigned v) const {
		ar & dir;
		ar & bw;
		boost::serialization::serialize(ar, *array_ptr, v);
	}

	template<typename Arc>
	void load(Arc& ar, const unsigned v) {
		isvirtual = true;
		array_ptr = std::make_shared<xarray<Ndim>>();
		ar & dir;
		ar & bw;
		boost::serialization::serialize(ar, *array_ptr, v);
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER();

	};

}

#endif /* XARRAY_HPP_ */
