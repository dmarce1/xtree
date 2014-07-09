/*
 * delayed_action.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: dmarce1
 */

#ifndef DELAYED_ACTION_HPP_
#define DELAYED_ACTION_HPP_

namespace xtree {

template<typename T>
class delayed_action {
	const T* input_reference;
	mutable T output;
	std::function<T(const T&)> action;
	mutable bool done;

	void do_action() const {
		if (!done) {
			output = action(*input_reference);
			done = true;
		}
	}
public:

	delayed_action() :
			input_reference(nullptr) {
		done = true;
	}
	template<class FuncType>
	delayed_action(const T* i, FuncType func) :
			input_reference(i), action(func) {
		done = false;
	}

	T get() {
		do_action();
		return std::move(output);
	}

	template<typename Arc>
	void load(Arc& ar, const unsigned v) {
		boost::serialization::serialize(ar, output, v);
		done = true;
	}

	template<typename Arc>
	void save(Arc& ar, const unsigned v) const {
		do_action();
		boost::serialization::serialize(ar, output, v);
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER();
};

}

#endif /* DELAYED_ACTION_HPP_ */
