/*
 * types.hpp
 *
 *  Created on: May 24, 2014
 *      Author: dmarce1
 */

#ifndef TYPES3_HPP_
#define TYPES3_HPP_

namespace xtree {

using iter_type = std::size_t;

template<int, int...>
struct int_seq;

template<int>
class silo_output;

template<typename >
struct int_seq_over2;

template<int, int, int>
class indexer;

template<int>
class location;

template<int>
class node_base;

template<typename T>
T factorial(T);

template<typename, int>
class node;

template<typename, int>
class tree;

template<int, int>
struct pow_;

template<int Ndim>
using child_index_type = indexer<Ndim, 2, 0>;

template<int Ndim>
using dir_type = indexer<Ndim, 3, -1>;

}

#endif /* TYPES_HPP_ */
