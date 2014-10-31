/*
 * exafmm.cpp
 *
 *  Created on: Oct 28, 2014
 *      Author: dmarce1
 */

#include "exafmm.hpp"
#include "tbb/cache_aligned_allocator.h"
#include <algorithm>

#define SIGN(m) ((m) < 0 ? -1.0 : 1.0)
//template class exafmm_kernel<5> ;
template class exafmm_kernel<5> ;

template<std::int64_t P>
std::valarray<std::valarray<complex>> exafmm_kernel<P>::M2L_interior(std::valarray<std::valarray<complex>> M, real dx,
		std::int64_t Nx, bool leaf, bool is_root) {
	const std::int64_t max_cnt = is_root ? Nx * Nx * Nx : 216;
	const std::int64_t min_dist = leaf ? 0 : 1;
	std::valarray<real> x(max_cnt), y(max_cnt), z(max_cnt);
	std::valarray<std::valarray<complex>> L(std::valarray<complex>(0.0, P * (P + 1) / 2), Nx * Nx * Nx);
	std::valarray<std::size_t> indexes(max_cnt);
	for (std::int64_t l0 = 0; l0 != Nx; ++l0) {
		const auto l1_min = is_root ? 0 : std::max(std::int64_t(0), ((l0 / 2) - 1) * 2);
		const auto l1_max = is_root ? Nx - 1 : std::min(Nx - 1, (((l0 / 2) + 1) * 2) + 1);
		for (std::int64_t k0 = 0; k0 != Nx; ++k0) {
			const auto k1_min = is_root ? 0 : std::max(std::int64_t(0), ((k0 / 2) - 1) * 2);
			const auto k1_max = is_root ? Nx - 1 : std::min(Nx - 1, (((k0 / 2) + 1) * 2) + 1);
			for (std::int64_t j0 = 0; j0 != Nx; ++j0) {
				const auto j1_min = is_root ? 0 : std::max(std::int64_t(0), ((j0 / 2) - 1) * 2);
				const auto j1_max = is_root ? Nx - 1 : std::min(Nx - 1, (((j0 / 2) + 1) * 2) + 1);

				std::int64_t cnt = 0;
				for (std::int64_t l1 = l1_min; l1 <= l1_max; ++l1) {
					for (std::int64_t k1 = k1_min; k1 <= k1_max; ++k1) {
						for (std::int64_t j1 = j1_min; j1 <= j1_max; ++j1) {
							std::int64_t max_dist = std::max(std::abs(k1 - k0), std::abs(l1 - l0));
							max_dist = std::max(std::abs(j1 - j0), max_dist);
							if (max_dist > min_dist) {
								indexes[cnt] = j1 + Nx * (k1 + Nx * l1);
								x[cnt] = dx * (j1 - j0);
								y[cnt] = dx * (k1 - k0);
								z[cnt] = dx * (l1 - l0);
								cnt++;
							}
						}
					}
				}
				auto this_L = M2L_vec(M[j0 + Nx * (k0 + Nx * l0)], x, y, z, cnt);
				for (std::int64_t n = 0; n != P * (P + 1) / 2; ++n) {
					for (std::int64_t i = 0; i != cnt; ++i) {
						L[indexes[i]][n] += this_L[n][i];
					}
				}
			}
		}
	}
	return L;
}


template<std::int64_t P>
void exafmm_kernel<P>::M2L(std::valarray<complex>& CiL, const std::valarray<complex> CjM, const std::valarray<real>& dist) {
	std::valarray<complex> Ynm(P * P);
	real rho, theta, phi;
	cart2sph(rho, theta, phi, dist);
	evalLocal(rho, theta, phi, Ynm);
	for (int j = 0; j != P; ++j) {
		for (int k = 0; k <= j; ++k) {
			int jk = j * j + j + k;
			int jks = j * (j + 1) / 2 + k;
			complex L = 0;
			for (int n = 0; n != P - j; ++n) {
				for (int m = -n; m < 0; ++m) {
					int nm = n * n + n + m;
					int nms = n * (n + 1) / 2 - m;
					int jknm = jk * P * P + nm;
					int jnkm = (j + n) * (j + n) + j + n + m - k;
					L += std::conj(CjM[nms]) * Cnm[jknm] * Ynm[jnkm];
				}
				for (int m = 0; m <= n; ++m) {
					int nm = n * n + n + m;
					int nms = n * (n + 1) / 2 + m;
					int jknm = jk * P * P + nm;
					int jnkm = (j + n) * (j + n) + j + n + m - k;
					L += CjM[nms] * Cnm[jknm] * Ynm[jnkm];
				}
			}
			CiL[jks] = L;
// printf( "%i %e %e\n", jks, L.real(), L.imag());
		}
	}
}

template<std::int64_t P>
void exafmm_kernel<P>::cart2sph(real& r, real& theta, real& phi, std::valarray<real> dist) {
	r = sqrt((dist * dist).sum()) * (1.0);
	// r = sqrt(x^2 + y^2 + z^2)
	if (r < EPS) {                                             // If r == 0
		theta = 0;
		//  theta can be anything so we set it to 0
	} else {                                                    // If r != 0
		theta = acos(dist[2] / r);
		//  theta = acos(z / r)
	}                                                   // End if for r == 0
	phi = atan2(dist[1], dist[0]);
}

template<std::int64_t P>
void exafmm_kernel<P>::M2M(std::valarray<complex>& CiM, const std::valarray<complex>& CjM,
		const std::valarray<real>& dist) {
	const complex I(0., 1.);
	std::valarray<complex> Ynm(P * P);

	real rho, theta, phi;
	cart2sph(rho, theta, phi, dist);
	evalMultipole(rho, theta, -phi, Ynm);
	for (std::int64_t j = 0; j != P; ++j) {
		for (std::int64_t k = 0; k <= j; ++k) {
			std::int64_t jk = j * j + j + k;
			std::int64_t jks = j * (j + 1) / 2 + k;
			complex M = 0;
			for (std::int64_t n = 0; n <= j; ++n) {
				for (std::int64_t m = -n; m <= std::min(k - 1, n); ++m) {
					if (j - n >= k - m) {
						std::int64_t jnkm = (j - n) * (j - n) + j - n + k - m;
						std::int64_t jnkms = (j - n) * (j - n + 1) / 2 + k - m;
						std::int64_t nm = n * n + n + m;
						M += CjM[jnkms] * std::pow(I, real(m - abs(m))) * Ynm[nm]
								* real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
					}
				}
				for (std::int64_t m = k; m <= n; ++m) {
					if (j - n >= m - k) {
						std::int64_t jnkm = (j - n) * (j - n) + j - n + k - m;
						std::int64_t jnkms = (j - n) * (j - n + 1) / 2 - k + m;
						std::int64_t nm = n * n + n + m;
						M += std::conj(CjM[jnkms]) * Ynm[nm] * real(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
					}
				}
			}
			CiM[jks] += M;
		}
	}
}

template<std::int64_t P>
void exafmm_kernel<P>::M2M_vec(std::valarray<std::valarray<complex>>& CiM,
		const std::valarray<std::valarray<complex>>& CjM, real x0, real y0, real z0, std::size_t vlen) {
	std::vector < std::vector < real >> Min_real(P * (P + 1) / 2, std::vector < real > (vlen));
	std::vector < std::vector < real >> Min_imag(P * (P + 1) / 2, std::vector < real > (vlen));
	std::vector < std::vector < real >> Mout_real(P * (P + 1) / 2, std::vector < real > (vlen));
	std::vector < std::vector < real >> Mout_imag(P * (P + 1) / 2, std::vector < real > (vlen));
	real rho, phi, theta;
	std::vector<real> Ynm_real(P * P);
	std::vector<real> Ynm_imag(P * P);

	for (std::size_t n = 0; n != P * (P + 1) / 2; ++n) {
		for (std::size_t i = 0; i != vlen; ++i) {
			Min_real[n][i] = CjM[i][n].real();
			Min_imag[n][i] = CjM[i][n].imag();
			Mout_real[n][i] = CiM[i][n].real();
			Mout_imag[n][i] = CiM[i][n].imag();
		}
	}
	rho = sqrt(x0 * x0 + y0 * y0 + z0 * z0);
	theta = acos(z0 / rho);
	phi = atan2(y0, x0);
	const complex I(0., 1.);                               // Imaginary unit
	phi = -phi;
	real x = std::cos(theta);                              // x = cos(theta)
	real y = std::sin(theta);                              // y = sin(theta)
	real fact = 1;                                   // Initialize 2 * m + 1
	real pn = 1;                        // Initialize Legendre polynomial Pn
	real rhom = 1;                                       // Initialize rho^m
	for (std::int64_t m = 0; m != P; ++m) {                     // Loop over m in Ynm
		complex eim = std::exp(I * real(m * phi));      //  exp(i * m * phi)
		real p = pn;                  //  Associated Legendre polynomial Pnm
		std::int64_t npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
		std::int64_t nmn = m * m;                          //  Index of Ynm for m < 0
		Ynm_real[npn] = rhom * p * prefactor[npn] * eim.real(); //  rho^m * Ynm for m > 0
		Ynm_imag[npn] = rhom * p * prefactor[npn] * eim.imag(); //  rho^m * Ynm for m > 0
		Ynm_real[nmn] = Ynm_real[npn]; //  Use conjugate relation for m < 0
		Ynm_imag[nmn] = -Ynm_imag[npn]; //  Use conjugate relation for m < 0
		real p1 = p;                                              //  Pnm-1
		p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
		rhom *= rho;                                              //  rho^m
		real rhon = rhom;                                         //  rho^n
		for (std::int64_t n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
			std::int64_t npm = n * n + n + m;             //   Index of Ynm for m > 0
			std::int64_t nmm = n * n + n - m;             //   Index of Ynm for m < 0
			Ynm_real[npm] = rhon * p * prefactor[npm] * eim.real();     //   rho^n * Ynm
			Ynm_imag[npm] = rhon * p * prefactor[npm] * eim.imag();     //   rho^n * Ynm
			Ynm_real[nmm] = Ynm_real[npm]; //   Use conjugate relation for m < 0
			Ynm_imag[nmm] = -Ynm_imag[npm]; //   Use conjugate relation for m < 0
			real p2 = p1;                                         //   Pnm-2
			p1 = p;                                               //   Pnm-1
			p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
			rhon *= rho;                                   //   Update rho^n
		}                                         //  End loop over n in Ynm
		pn = -pn * fact * y;                                      //  Pn
		fact += 2;                                             //  2 * m + 1
	}                                              // End loop over m in Ynm

	for (std::int64_t j = 0; j != P; ++j) {
		for (std::int64_t k = 0; k <= j; ++k) {
			std::int64_t jk = j * j + j + k;
			std::int64_t jks = j * (j + 1) / 2 + k;
			for (std::int64_t n = 0; n <= j; ++n) {
				for (std::int64_t m = -n; m <= n; ++m) {
					if (j - n >= std::abs(k - m)) {
#pragma simd
#pragma vector aligned
						for (std::int64_t i = 0; i != vlen; ++i) {
							const std::int64_t jnkm = (j - n) * (j - n) + j - n + k - m;
							const std::int64_t jnkms = (j - n) * (j - n + 1) / 2 + std::abs(k - m);
							const std::int64_t nm = n * n + n + m;
							real tmp = real(ODDEVEN(n+(m>=k?m+k:0)) * Anm[nm] * Anm[jnkm] / Anm[jk]);
							real tmp_r = Min_real[jnkms][i] * Ynm_real[nm]
									- Min_imag[jnkms][i] * SIGN(k - m) * Ynm_imag[nm];
							real tmp_i = Min_imag[jnkms][i] * SIGN( k - m) * Ynm_real[nm]
									+ Min_real[jnkms][i] * Ynm_imag[nm];
							tmp *= (m >= 0) || (-m % 2 == 0) ? 1.0 : -1.0;
							tmp_r *= tmp;
							tmp_i *= tmp;
							Mout_real[jks][i] += tmp_r;
							Mout_imag[jks][i] += tmp_i;
						}
					}
				}
			}
		}
	}
	for (std::size_t i = 0; i != vlen; ++i) {
		for (std::size_t n = 0; n != P * (P + 1) / 2; ++n) {
			CiM[i][n].real(Mout_real[n][i]);
			CiM[i][n].imag(Mout_imag[n][i]);
		}
	}
}

template<std::int64_t P>
std::valarray<std::valarray<complex>> exafmm_kernel<P>::M2L_vec(const std::valarray<complex> CjM,
		const std::valarray<real>& x0, const std::valarray<real>& y0, const std::valarray<real>& z0, std::size_t vlen) {
	std::vector<real> rho(vlen), theta(vlen), phi(vlen);
	std::vector < std::vector < real >> Ynm_real(P * P, std::vector < real > (vlen));
	std::vector < std::vector < real >> Ynm_imag(P * P, std::vector < real > (vlen));
	std::vector<real> L_real(vlen);
	std::vector<real> L_imag(vlen);
	std::valarray<std::valarray<complex>> CiL(std::valarray<complex>(vlen), P * (P + 1) / 2);
#pragma simd
#pragma vector aligned
	for (std::size_t i = 0; i != vlen; ++i) {
		rho[i] = sqrt(x0[i] * x0[i] + y0[i] * y0[i] + z0[i] * z0[i]);
		theta[i] = acos(z0[i] / rho[i]);
		phi[i] = atan2(y0[i], x0[i]);
	}
#pragma simd
#pragma vector aligned
	for (std::size_t i = 0; i != vlen; ++i) {
		const complex I(0., 1.);                               // Imaginary unit
		real x = std::cos(theta[i]);                           // x = cos(theta)
		real y = std::sin(theta[i]);                           // y = sin(theta)
		real fact = 1;                                   // Initialize 2 * m + 1
		real pn = 1;                        // Initialize Legendre polynomial Pn
		real rhom = 1.0 / rho[i];                       // Initialize rho^(-m-1)
#pragma novector                                  //  rho^(-n-1)
		for (std::int64_t m = 0; m != P; ++m) {                     // Loop over m in Ynm
			real eim_real = cos(m * phi[i]);      //  exp(i * m * phi)
			real eim_imag = sin(m * phi[i]);      //  exp(i * m * phi)
			real p = pn;                  //  Associated Legendre polynomial Pnm
			std::int64_t npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
			std::int64_t nmn = m * m;                          //  Index of Ynm for m < 0
			Ynm_real[npn][i] = rhom * p * prefactor[npn] * eim_real; //  rho^(-m-1) * Ynm for m > 0
			Ynm_imag[npn][i] = rhom * p * prefactor[npn] * eim_imag; //  rho^(-m-1) * Ynm for m > 0
			Ynm_real[nmn][i] = Ynm_real[npn][i];
			Ynm_imag[nmn][i] = -Ynm_imag[npn][i];
			real p1 = p;                                              //  Pnm-1
			p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
			rhom /= rho[i];                                       //  rho^(-m-1)
			real rhon = rhom;
#pragma novector                                  //  rho^(-n-1)
			for (std::int64_t n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
				std::int64_t npm = n * n + n + m;             //   Index of Ynm for m > 0
				std::int64_t nmm = n * n + n - m;             //   Index of Ynm for m < 0
				Ynm_real[npm][i] = rhon * p * prefactor[npm] * eim_real; //   rho^n * Ynm for m > 0
				Ynm_imag[npm][i] = rhon * p * prefactor[npm] * eim_imag; //   rho^n * Ynm for m > 0
				Ynm_real[nmm][i] = Ynm_real[npm][i];
				Ynm_imag[nmm][i] = -Ynm_imag[npm][i];
				real p2 = p1;                                         //   Pnm-2
				p1 = p;                                               //   Pnm-1
				p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
				rhon /= rho[i];                                  //   rho^(-n-1)
			}                                         //  End loop over n in Ynm
			pn = -pn * fact * y;                                      //  Pn
			fact += 2.0;                                           //  2 * m + 1
		}                                              // End loop over m in Ynm
	}
	for (std::int64_t j = 0; j != P; ++j) {
		for (std::int64_t k = 0; k <= j; ++k) {
			const std::int64_t jk = j * j + j + k;
			const std::int64_t jks = j * (j + 1) / 2 + k;
#pragma simd
#pragma vector aligned
			for (std::size_t i = 0; i != vlen; ++i) {
				L_real[i] = 0.0;
				L_imag[i] = 0.0;
			}
			for (std::int64_t n = 0; n != P - j; ++n) {
				for (std::int64_t m = -n; m <= n; ++m) {
					const std::int64_t nm = n * n + n + m;
					const std::int64_t nms = n * (n + 1) / 2 + std::abs(m);
					const std::int64_t jknm = jk * P * P + nm;
					const std::int64_t jnkm = (j + n) * (j + n) + j + n + m - k;
					const real Cnm_r = Cnm_real[jknm];
					const real Cnm_i = Cnm_imag[jknm];
					const real M_r = CjM[nms].real();
					const real M_i = CjM[nms].imag() * SIGN(m);
#pragma simd
#pragma vector aligned
					for (std::size_t i = 0; i != vlen; ++i) {
						real tmp_real = Cnm_r * Ynm_real[jnkm][i] - Cnm_i * Ynm_imag[jnkm][i];
						real tmp_imag = Cnm_r * Ynm_imag[jnkm][i] + Cnm_i * Ynm_real[jnkm][i];
						L_real[i] += tmp_real * M_r - tmp_imag * M_i;
						L_imag[i] += tmp_real * M_i + tmp_imag * M_r;
					}
				}
			}
#pragma simd
			for (std::size_t i = 0; i != vlen; ++i) {
				CiL[jks][i].real(L_real[i]);
				CiL[jks][i].imag(L_imag[i]);
			}
		}
	}
	return std::move(CiL);
}

template<std::int64_t P>
void exafmm_kernel<P>::L2L(std::valarray<complex>& CiL, const std::valarray<complex>& CjL,
		const std::valarray<real>& dist) {
	const complex I(0., 1.);
	std::valarray<complex> Ynm(P * P);
	real rho, theta, phi;
	cart2sph(rho, theta, phi, dist);
	evalMultipole(rho, theta, phi, Ynm);
	for (std::int64_t j = 0; j != P; ++j) {
		for (std::int64_t k = 0; k <= j; ++k) {
			std::int64_t jk = j * j + j + k;
			std::int64_t jks = j * (j + 1) / 2 + k;
			complex L = 0;
			for (std::int64_t n = j; n != P; ++n) {

				for (std::int64_t m = j + k - n; m < 0; ++m) {
					std::int64_t jnkm = (n - j) * (n - j) + n - j + m - k;
					std::int64_t nm = n * n + n - m;
					std::int64_t nms = n * (n + 1) / 2 - m;
					L += std::conj(CjL[nms]) * Ynm[jnkm] * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
				}
				for (std::int64_t m = 0; m <= n; ++m) {
					if (n - j >= abs(m - k)) {
						std::int64_t jnkm = (n - j) * (n - j) + n - j + m - k;
						std::int64_t nm = n * n + n + m;
						std::int64_t nms = n * (n + 1) / 2 + m;
						L += CjL[nms] * std::pow(I, real(m - k - abs(m - k))) * Ynm[jnkm] * Anm[jnkm] * Anm[jk]
								/ Anm[nm];
					}
				}
			}
			CiL[jks] = L;
		}
	}
}

template<std::int64_t P>
void exafmm_kernel<P>::evalMultipole(real rho, real theta, real phi, std::valarray<complex>& Ynm) {
	const complex I(0., 1.);                               // Imaginary unit
	real x = std::cos(theta);                              // x = cos(theta)
	real y = std::sin(theta);                              // y = sin(theta)
	real fact = 1;                                   // Initialize 2 * m + 1
	real pn = 1;                        // Initialize Legendre polynomial Pn
	real rhom = 1;                                       // Initialize rho^m
	for (std::int64_t m = 0; m != P; ++m) {                     // Loop over m in Ynm
		complex eim = std::exp(I * real(m * phi));      //  exp(i * m * phi)
		real p = pn;                  //  Associated Legendre polynomial Pnm
		std::int64_t npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
		std::int64_t nmn = m * m;                          //  Index of Ynm for m < 0
		Ynm[npn] = rhom * p * prefactor[npn] * eim; //  rho^m * Ynm for m > 0
		Ynm[nmn] = std::conj(Ynm[npn]); //  Use conjugate relation for m < 0
		real p1 = p;                                              //  Pnm-1
		p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
		rhom *= rho;                                              //  rho^m
		real rhon = rhom;                                         //  rho^n
		for (std::int64_t n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
			std::int64_t npm = n * n + n + m;             //   Index of Ynm for m > 0
			std::int64_t nmm = n * n + n - m;             //   Index of Ynm for m < 0
			Ynm[npm] = rhon * p * prefactor[npm] * eim;     //   rho^n * Ynm
			Ynm[nmm] = std::conj(Ynm[npm]); //   Use conjugate relation for m < 0
			real p2 = p1;                                         //   Pnm-2
			p1 = p;                                               //   Pnm-1
			p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
			rhon *= rho;                                   //   Update rho^n
		}                                         //  End loop over n in Ynm
		pn = -pn * fact * y;                                      //  Pn
		fact += 2;                                             //  2 * m + 1
	}                                              // End loop over m in Ynm
}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
template<std::int64_t P>
void exafmm_kernel<P>::evalLocal(real rho, real theta, real phi, std::valarray<complex>& Ynm) {
	const complex I(0., 1.);                               // Imaginary unit
	real x = std::cos(theta);                              // x = cos(theta)
	real y = std::sin(theta);                              // y = sin(theta)
	real fact = 1;                                   // Initialize 2 * m + 1
	real pn = 1;                        // Initialize Legendre polynomial Pn
	real rhom = 1.0 / rho;                          // Initialize rho^(-m-1)
	for (std::int64_t m = 0; m != P; ++m) {                     // Loop over m in Ynm
		complex eim = std::exp(I * real(m * phi));      //  exp(i * m * phi)
		real p = pn;                  //  Associated Legendre polynomial Pnm
		std::int64_t npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
		std::int64_t nmn = m * m;                          //  Index of Ynm for m < 0
		Ynm[npn] = rhom * p * prefactor[npn] * eim; //  rho^(-m-1) * Ynm for m > 0
		Ynm[nmn] = std::conj(Ynm[npn]); //  Use conjugate relation for m < 0
		real p1 = p;                                              //  Pnm-1
		p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
		rhom /= rho;                                          //  rho^(-m-1)
		real rhon = rhom;                                     //  rho^(-n-1)
		for (std::int64_t n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
			std::int64_t npm = n * n + n + m;             //   Index of Ynm for m > 0
			std::int64_t nmm = n * n + n - m;             //   Index of Ynm for m < 0
			Ynm[npm] = rhon * p * prefactor[npm] * eim; //   rho^n * Ynm for m > 0
			Ynm[nmm] = std::conj(Ynm[npm]); //   Use conjugate relation for m < 0
			real p2 = p1;                                         //   Pnm-2
			p1 = p;                                               //   Pnm-1
			p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
			rhon /= rho;                                     //   rho^(-n-1)
		}                                         //  End loop over n in Ynm
		pn = -pn * fact * y;                                      //  Pn
		fact += 2;                                             //  2 * m + 1
	}                                              // End loop over m in Ynm
}

template<std::int64_t P>
exafmm_kernel<P>::exafmm_kernel() {
	const complex I(0., 1.);                               // Imaginary unit

	factorial[0] = 1;                                // Initialize factorial
	for (std::int64_t n = 1; n != P; ++n) {                              // Loop to P
		factorial[n] = factorial[n - 1] * n;                        //  n!
	}                                                       // End loop to P

	for (std::int64_t n = 0; n != P; ++n) {                     // Loop over n in Anm
		for (std::int64_t m = -n; m <= n; ++m) {               //  Loop over m in Anm
			std::int64_t nm = n * n + n + m;                        //   Index of Anm
			std::int64_t nabsm = abs(m);                                     //   |m|
			real fnmm = 1.0;                        //   Initialize (n - m)!
			for (std::int64_t i = 1; i <= n - m; ++i)
				fnmm *= i;                  //   (n - m)!
			real fnpm = 1.0;                        //   Initialize (n + m)!
			for (std::int64_t i = 1; i <= n + m; ++i)
				fnpm *= i;                  //   (n + m)!
			real fnma = 1.0;                      //   Initialize (n - |m|)!
			for (std::int64_t i = 1; i <= n - nabsm; ++i)
				fnma *= i;              //   (n - |m|)!
			real fnpa = 1.0;                      //   Initialize (n + |m|)!
			for (std::int64_t i = 1; i <= n + nabsm; ++i)
				fnpa *= i;              //   (n + |m|)!
			prefactor[nm] = std::sqrt(fnma / fnpa); //   sqrt( (n - |m|)! / (n + |m|)! )
			Anm[nm] = ODDEVEN(n) / std::sqrt(fnmm * fnpm); //   (-1)^n / sqrt( (n + m)! / (n - m)! )
		}                                         //  End loop over m in Anm
	}                                              // End loop over n in Anm

	for (std::int64_t j = 0, jk = 0, jknm = 0; j != P; ++j) { // Loop over j in Cjknm
		for (std::int64_t k = -j; k <= j; ++k, ++jk) {       //  Loop over k in Cjknm
			for (std::int64_t n = 0, nm = 0; n != P; ++n) { //   Loop over n in Cjknm
				for (std::int64_t m = -n; m <= n; ++m, ++nm, ++jknm) { //    Loop over m in Cjknm
					if (j + n < P) {
						const std::int64_t jnkm = (j + n) * (j + n) + j + n + m - k; //     Index C_{j+n}^{m-k}
						Cnm[jknm] = std::pow(I, real(abs(k - m) - abs(k) - abs(m))) //     Cjknm
						* real(ODDEVEN(j) * Anm[nm] * Anm[jk] / Anm[jnkm]);
						Cnm_real[jknm] = Cnm[jknm].real();
						Cnm_imag[jknm] = Cnm[jknm].imag();
					}                         //    End loop over m in Cjknm
				}
			}             //   End loop over n in Cjknm
		}                                    //  End loop over in k in Cjknm
	}                                         // End loop over in j in Cjknm
}

template<std::int64_t P>
std::valarray<real> exafmm_kernel<P>::factorial(P);

template<std::int64_t P>
std::valarray<real> exafmm_kernel<P>::prefactor(P * P);

template<std::int64_t P>
std::valarray<real> exafmm_kernel<P>::Anm(P * P);

template<std::int64_t P>
std::valarray<complex> exafmm_kernel<P>::Cnm(P * P * P * P);

template<std::int64_t P>
std::valarray<real> exafmm_kernel<P>::Cnm_real(P * P * P * P);

template<std::int64_t P>
std::valarray<real> exafmm_kernel<P>::Cnm_imag(P * P * P * P);

