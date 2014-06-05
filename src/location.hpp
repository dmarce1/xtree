/*
 * location.hpp
 *
 *  Created on: May 22, 2014
 *      Author: dmarce1
 */

#ifndef LOCATION_HPP_
#define LOCATION_HPP_

namespace xtree {

template<int Ndim>
class location {
private:
	int level;
	std::array<int, Ndim> loc;
public:
	location() {
		level = 0;
		std::fill(loc.begin(), loc.end(), 0);
	}
	virtual ~location() {
	}
	location(const location& l) {
		level = l.level;
		loc = l.loc;
	}
	location& operator=(const location& l) {
		level = l.level;
		loc = l.loc;
		return *this;
	}
	template<typename Arc>
	void serialize(Arc& ar, const int v) {
		ar & level;
		ar & loc;
	}
	dir_type<Ndim> relative_direction_to(const location<Ndim>& loc2) {
		assert(loc2.level == level);
		dir_type<Ndim> dir;
		for (int d = 0; d < Ndim; d++) {
			const int dif = loc2.loc[d] - loc[d];
			if (dif > 0) {
				dir[d] = +1;
			} else if (dif < 0) {
				dir[d] = -1;
			} else {
				dir[d] = 0;
			}
		}
		return dir;
	}
	location get_neighbor(const indexer<Ndim, 3, -1>& dir) {
		location rloc;
		for (int i = 0; i < Ndim; i++) {
			rloc.loc[i] = loc[i] + dir[i];
		}
		rloc.level = level;
		return rloc;
	}
	int get_level() const {
		return level;
	}
	void set_level(int i) {
		level = i;
	}
	std::array<int, Ndim> get_location() const {
		return loc;
	}
	int get_location(int i) const {
		return loc[i];
	}
	location& shift(int d, int amount) {
		loc[d] += amount;
		return *this;
	}
	void set_location(const std::array<int, Ndim>& p) {
		loc = p;
	}
	void set_location(int i, int j) {
		loc[i] = j;
	}
	std::array<double, Ndim> get_position() const {
		std::array<double, Ndim> pos;
		const double dx = get_dx();
		for (int i = 0; i < Ndim; i++) {
			pos[i] = loc[i] * dx;
		}
		return pos;
	}
	double get_dx() const {
		return 1.0 / double(1 << level);
	}
	location get_parent() const {
		location p = *this;
		p.level--;
		p.loc = p.loc / 2;
		return p;
	}
	indexer<Ndim, 2> this_child_index() const {
		indexer<Ndim, 2> ci;
		for (int i = 0; i < Ndim; i++) {
			ci[i] = loc[i] & 1;
		}
		return ci;
	}
	location get_child(const indexer<Ndim, 2>& ci) const {
		location c = *this;
		c.level++;
		c.loc = c.loc * 2;
		for (int i = 0; i < Ndim; i++) {
			c.loc[i] += ci[i];
		}
		return c;
	}
};
}
#endif /* LOCATION_HPP_ */
