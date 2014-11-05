/*
 * cube_poles.hpp
 *
 *  Created on: Jun 28, 2014
 *      Author: dmarce1
 */

#ifndef CUBE_POLES_HPP_
#define CUBE_POLES_HPP_

template<std::int64_t P>
class cube_poles {
	static constexpr int N = 100;
	std::valarray<real> M;
	std::valarray<real> NN;
public:
	std::valarray<real> get_M(double dx) {
		std::valarray < real > rc(M.size());
		for (auto i = 0; i < M.size(); i++) {
			rc[i] = pow(dx, NN[i]) * M[i];
		}
		return rc;
	}
	cube_poles() {
		exafmm_kernel < P > exafmm;
		std::valarray < real > Y;
		std::array < real, N + 1 > w;
		Y.resize(P * P);
		M.resize(P * P);
		const real dx = 1.0 / real(N);
		const real dx3 = dx * dx * dx;
		Y = 0.0;
		real r, theta, phi;
		for (int i = 0; i <= N; i++) {
			if (i == 0 || i == N) {
				w[i] == real(1.0 / 3.0);
			} else if (i % 2 == 1) {
				w[i] = real(2.0 / 3.0);
			} else {
				w[i] = real(4.0 / 3.0);
			}
		}
		for (int i = 0; i <= N; i++) {
			for (int j = 0; j <= N; j++) {
				for (int k = 0; k <= N; k++) {
					std::valarray < real > X(3);
					X[0] = (real(i) - real(N) / 2.0) * dx * 2.0;
					X[1] = (real(k) - real(N) / 2.0) * dx * 2.0;
					X[2] = (real(j) - real(N) / 2.0) * dx * 2.0;
					std::valarray < real > Ynm(P * P);
					std::valarray < real > Ynm_theta(P * P);
					exafmm.cart2sph(r, theta, phi, X);
					if (r > 1.0) {
						exafmm.evalMultipole(r, theta, phi, Ynm);
						Y += w[i] * w[j] * w[k] * Ynm;
					}
				}
			}
		}
		Y *= dx3;
		Y[0] = 1.0;
		M = Y;
		NN.resize(N * N);
		for (int n = 0; n != P; ++n) {
			for (int m = 0; m <= n; ++m) {
				NN[n * n + n + m] = NN[n * n + n - m] = real(n + 3);
			}
		}
	}
};

#endif /* CUBE_POLES_HPP_ */
