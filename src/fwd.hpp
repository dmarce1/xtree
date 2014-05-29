/*
 * types.hpp
 *
 *  Created on: May 24, 2014
 *      Author: dmarce1
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_

namespace xtree {

template<int, int...>
struct int_seq;

template<int, int, int>
class indexer;

template<int>
class location;

template<typename, int>
class node;

template<typename, int>
class tree;

template<int, int>
struct pow_;

template<typename, int...>
class vector;

template<int>
class array_index;

template<typename, int>
class grid_base;

template<typename, typename>
class subgrid;

template<typename, typename, int>
class grid;

template<typename, typename, int>
class bgrid;

template<int>
class grid_index;

template<typename T>
bool if_boolean_expression(T);

}

#endif /* TYPES_HPP_ */
