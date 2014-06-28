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
	std::valarray<complex> M;
	std::valarray<complex> NN;
public:
	std::valarray<complex> get_M(double dx) {
		return M * std::pow(complex(dx), NN);
	}
	cube_poles() {
		exafmm_kernel<P> exafmm;
		std::valarray<complex> Y;
		std::array<complex, N + 1> w;
		Y.resize(P * P);
		M.resize(P * (P + 1) / 2);
		const real dx = 1.0 / real(N);
		const real dx3 = dx * dx * dx;
		Y = 0.0;
		real r, theta, phi;
		for (int i = 0; i <= N; i++) {
			if (i == 0 || i == N) {
				w[i] == complex(1.0 / 3.0, 0.0);
			} else if (i % 2 == 1) {
				w[i] = complex(2.0 / 3.0, 0.0);
			} else {
				w[i] = complex(4.0 / 3.0, 0.0);
			}
		}
		for (int i = 0; i <= N; i++) {
			for (int j = 0; j <= N; j++) {
				for (int k = 0; k <= N; k++) {
					std::valarray<real> X(3);
					X[0] = (real(i) - real(N) / 2.0) * dx * 2.0;
					X[1] = (real(k) - real(N) / 2.0) * dx * 2.0;
					X[2] = (real(j) - real(N) / 2.0) * dx * 2.0;
					std::valarray<complex> Ynm(P * P);
					std::valarray<complex> Ynm_theta(P * P);
					exafmm.cart2sph(r, theta, phi, X);
					if (r > 1.0) {
						exafmm.evalMultipole(r, theta, phi, std::begin(Ynm), std::begin(Ynm_theta));
						Y += w[i] * w[j] * w[k] * Ynm;
					}
				}
			}
		}
		Y[0] = 1.0;
		Y *= dx3;
		for (int n = 0; n != P; ++n) {
			for (int m = 0; m <= n; ++m) {
				M[n * (n + 1) / 2 + m] = 0.5 * (Y[n * n + n + m] + Y[n * n + n - m]);
			}
		}
		NN.resize(N * (N + 1) / 2);
		for (int n = 0; n != P; ++n) {
			for (int m = 0; m <= n; ++m) {
				NN[n * (n + 1) / 2 + m] = double(n);
				M[n * (n + 1) / 2 + m] = 0.5 * (Y[n * n + n + m] + Y[n * n + n - m]);
			}
		}
	}
};

#endif /* CUBE_POLES_HPP_ */
