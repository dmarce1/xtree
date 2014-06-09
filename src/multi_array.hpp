/*
 * multi_array.hpp
 *
 *  Created on: Jun 9, 2014
 *      Author: dmarce1
 */

#ifndef MULTI_ARRAY_HPP_
#define MULTI_ARRAY_HPP_

#include <boost/multi_array.hpp>

namespace boost {
namespace serialization {
template<class Archive, class T, std::size_t Ndim>
void save(Archive & ar, const multi_array<T, Ndim>& t, const unsigned int) {
}
template<class Archive, class T, std::size_t Ndim>
void load(Archive & ar, multi_array<T, Ndim>& t, const unsigned int) {
}
template<class Archive, class T, std::size_t Ndim>
void inline serialize(Archive & ar, multi_array<T, Ndim>& t, const unsigned int file_version) {
	split_free(ar, t, file_version);
}

}
}

#endif /* MULTI_ARRAY_HPP_ */
