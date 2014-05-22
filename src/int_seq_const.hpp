/*
 * int_seq_const.hpp
 *
 *  Created on: May 21, 2014
 *      Author: dmarce1
 */

#ifndef INT_SEQ_CONST_HPP_
#define INT_SEQ_CONST_HPP_

namespace xtree {

template<int N, int Ndim>
struct int_seq_const {
	static constexpr int dim() {
		return Ndim;
	}
	static constexpr int get(int i) {
		return N;
	}
};


}




#endif /* INT_SEQ_CONST_HPP_ */
