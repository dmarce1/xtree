/*
 * types.hpp
 *
 *  Created on: May 24, 2014
 *      Author: dmarce1
 */

#ifndef TYPES3_HPP_
#define TYPES3_HPP_

namespace xtree {

template<int, int...>
struct int_seq;

template<int, int>
class silo_output;

template<typename >
struct int_seq_over2;

template<int, int, int>
class indexer;

template<int>
class location;

template<int>
class node_base;

template<typename, int>
class node;

template<typename, int>
class tree;

template<int, int>
struct pow_;


template<typename T>
bool if_boolean_expression(T);

}

#endif /* TYPES_HPP_ */
