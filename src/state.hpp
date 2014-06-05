/*
 * state.hpp
 *
 *  Created on: Jun 5, 2014
 *      Author: dmarce1
 */

#ifndef STATE_HPP_
#define STATE_HPP_

namespace xtree {

template<int Ndim, int Nf = Ndim + 2>
class state: public std::array<double, Nf> {
public:
	static constexpr double gamma = 5.0 / 3.0;
public:
	state() = default;
	state(const state&) = default;
	state(state&&) = default;
	state& operator=(const state&) = default;
	state& operator=(state&&) = default;

	template<typename T, std::size_t ...Size, template<typename, std::size_t...> class Container >
	state( const Container<T,Size...>& u ) {
		std::copy( u.begin(), u.end(), std::array<double,Nf>::begin());
	}

	void set_momentum(int d, double m) {
		std::array<double,Nf>::operator[](d) = m;
	}

	void set_density( double d) {
		std::array<double,Nf>::operator[](Ndim) = d;
	}

	void set_energy(double e) {
		std::array<double,Nf>::operator[](Ndim+1) = e;
	}

	double get_momentum(int d) const {
		return std::array<double,Nf>::operator[](d);
	}

	std::array<double, Ndim> get_momentum() const {
		std::array<double, Ndim> m;
		for( int i = 0; i < Ndim; i++) {
			m[i] = get_momentum(i);
		}
		return m;
	}

	double get_density() const {
		return std::array<double,Nf>::operator[](Ndim);
	}

	double get_energy() const {
		return std::array<double,Nf>::operator[](Ndim+1);
	}

	double get_velocity(int d) const {
		return get_momentum(d) / get_density();
	}

	std::array<double,Ndim> get_velocity() const {
		return get_momentum() / get_density();

	}

	double get_kinetic_energy() const {
		const std::array<double,Ndim> mom = get_momentum();
		double m2 = std::inner_product(mom.begin(), mom.end(), mom.begin(),0.0);
		return 0.5 * m2 / get_density();
	}

	double get_internal_energy() const {
		return std::max(get_energy() - get_kinetic_energy(), 0.0);
	}

	double get_pressure() const {
		return (gamma - 1.0) * get_internal_energy();
	}

	double get_sound_speed() const {
		return sqrt( gamma * get_pressure() / get_density() );
	}

	double get_signal_speed(int d ) const {
		return get_sound_speed() + std::abs(get_velocity(d));
	}

};

}

#endif /* STATE_HPP_ */
