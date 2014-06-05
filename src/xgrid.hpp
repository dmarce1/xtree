/*
 * xgrid.hpp
 *
 *  Created on: May 29, 2014
 *      Author: dmarce1
 */

#ifndef SUBGRID43_HPP_
#define SUBGRID43_HPP_

#include "fwd.hpp"

namespace xtree {

template<typename T, typename Dims>
class xgrid: public grid_base<T, Dims::dim()> {
public:
	static constexpr int Ndim = Dims::dim();
	static constexpr int Size = Dims::size() / (1 << Ndim);
	using base_type = grid_base<T, Ndim>;
	using index_type = typename grid_base<T, Ndim>::index_type;
	using grid_type = grid_base<T,Ndim>;
private:
	const child_index_type<Ndim> xtant;
	const vector<int, Ndim> origin;
	std::shared_ptr<const grid_type> parent;
	std::shared_ptr<vector<T, Size>> data_ptr;
public:
	xgrid() {
	}
	xgrid(std::shared_ptr<const grid_type>& _parent, const child_index_type<Ndim>& _xtant) :
			parent(_parent), xtant(_xtant), origin((xtant.to_vector() * Dims::to_vector()) / 2) {
	}
	virtual ~xgrid() {
	}
	xgrid(xgrid<T, Dims> && x) :
			xtant(std::move(x.xtant)), origin(std::move(x.origin)) {
		parent = std::move(x.parent);
		data_ptr = std::move(x.data_ptr);
	}
	template<typename Arc>
	void load(Arc& ar, const int v) {
		data_ptr = std::make_shared<vector<T, Size>>();
		data_ptr->load(ar, v);
	}
	template<typename Arc>
	void save(Arc& ar, const int v) const {
		if (parent) {
			for (grid_index<Ndim> ci(Dims::to_vector() - 1); !ci.end(); ++ci) {
				ar & (*this)[ci];
			}
		} else {
			data_ptr->save(ar, v);
		}

	}
	const T& get(const index_type& i) const {
		if (parent) {
			return parent->data[base_type::template vector_to_index<Dims>(i + origin)];
		} else {

		}
	}
};

}

#endif /* SUBGRID_HPP_ */
