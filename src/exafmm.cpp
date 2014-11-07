/*
 * exafmm.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: dmarce1
 */

#define EXAFMM_CPP
#include "exafmm.hpp"
#include "offload.h"
#include <array>

#define ODDEVEN(n) real((((n) & 1) == 1) ? -1 : 1)

#define COMPLEX_MULT( ar, ai, br, bi, cr, ci) \
		(ar) = (br)*(cr)-(bi)*(ci);           \
		(ai) = (br)*(ci)+(bi)*(cr)

#define COMPLEX_MULT_ADD( ar, ai, br, bi, cr, ci) \
		(ar) += (br)*(cr)-(bi)*(ci);           \
		(ai) += (br)*(ci)+(bi)*(cr)

#define SGN(i) real((i) > 0 ? 1 : ((i)<0 ? -1 : 0))

template class exafmm_kernel<1> ;
template class exafmm_kernel<2> ;
template class exafmm_kernel<3> ;
template class exafmm_kernel<4> ;
template class exafmm_kernel<5> ;
template class exafmm_kernel<6> ;
template class exafmm_kernel<7> ;
template class exafmm_kernel<8> ;
template class exafmm_kernel<9> ;
template class exafmm_kernel<10> ;

namespace mic {
template<class T>
__attribute__ (( target (mic)))
inline T max(const T& a, const T& b) {
	return a >= b ? a : b;
}
template<class T>
__attribute__ (( target (mic)))
inline T min(const T& a, const T& b) {
	return a <= b ? a : b;
}
}

template<std::int64_t P>
std::valarray<std::valarray<real>> exafmm_kernel<P>::M2L_interior(const std::valarray<std::valarray<real>>& M, real dx,
		std::int64_t Nx, bool leaf, bool is_root) {
	constexpr auto P2 = P * P;
	const std::size_t sz = Nx * Nx * Nx;
	const std::size_t buffer_size = P2 * sz;
	std::valarray<std::valarray<real>> L(std::valarray<real>(sz), P2);
	real* m_buffer = new real[buffer_size];
	real* l_buffer = new real[buffer_size];

	std::size_t cnt = 0;
	for (std::size_t p = 0; p != P2; ++p) {
		for (std::size_t i = 0; i != sz; ++i) {
			m_buffer[cnt++] = M[p][i];
		}
	}
	const real *_prefactor = prefactor;
	const real *_Cnm_r = Cnm_r;
	const real *_Cnm_i = Cnm_i;
#pragma offload target(mic) in(m_buffer : length(buffer_size)) in(_prefactor : length(P*P)) in(_Cnm_r : length(P*P*P*P)) in(_Cnm_i : length(P*P*P*P)) in(sz, is_root, leaf, Nx) out(l_buffer: length(buffer_size))
	{
		for (std::size_t i = 0; i != P2 * sz; ++i) {
			l_buffer[i] = real(0.0);
		}
		std::int64_t max_cnt = is_root ? sz : 216;
		if ((max_cnt * sizeof(real)) % 64 != 0) {
			const auto d = 64 / sizeof(real);
			max_cnt = ((max_cnt / d) + 1) * d;
		}
		const std::int64_t min_dist = leaf ? 0 : 1;
		auto* buffer = new real[(3 + 1 + 1 + P2) * max_cnt];
		auto* indexes = new std::size_t[max_cnt];
		auto* X = buffer + 0;
		auto* L_r = buffer + 3 * max_cnt;
		auto* L_i = buffer + 4 * max_cnt;
		auto* Ynm = buffer + 5 * max_cnt;
		for (std::int64_t l0 = 0; l0 != Nx; ++l0) {
			const auto l1_min = is_root ? 0 : mic::max(std::int64_t(0), ((l0 / 2) - 1) * 2);
			const auto l1_max = is_root ? Nx - 1 : mic::min(Nx - 1, (((l0 / 2) + 1) * 2) + 1);
			for (std::int64_t k0 = 0; k0 != Nx; ++k0) {
				const auto k1_min = is_root ? 0 : mic::max(std::int64_t(0), ((k0 / 2) - 1) * 2);
				const auto k1_max = is_root ? Nx - 1 : mic::min(Nx - 1, (((k0 / 2) + 1) * 2) + 1);
				for (std::int64_t j0 = 0; j0 != Nx; ++j0) {
					const auto j1_min = is_root ? 0 : mic::max(std::int64_t(0), ((j0 / 2) - 1) * 2);
					const auto j1_max = is_root ? Nx - 1 : mic::min(Nx - 1, (((j0 / 2) + 1) * 2) + 1);
					std::int64_t cnt = 0;
					for (std::int64_t l1 = l1_min; l1 <= l1_max; ++l1) {
						for (std::int64_t k1 = k1_min; k1 <= k1_max; ++k1) {
							for (std::int64_t j1 = j1_min; j1 <= j1_max; ++j1) {
								std::int64_t max_dist = mic::max(std::abs(k1 - k0), std::abs(l1 - l0));
								max_dist = mic::max(std::abs(j1 - j0), max_dist);
								if (max_dist > min_dist) {
									indexes[cnt] = j1 + Nx * (k1 + Nx * l1);
									X[0 * max_cnt + cnt] = dx * (j1 - j0);
									X[1 * max_cnt + cnt] = dx * (k1 - k0);
									X[2 * max_cnt + cnt] = dx * (l1 - l0);
									cnt++;
								}
							}
						}
					}
#pragma vector aligned
#pragma simd
					for (std::size_t i = 0; i != cnt; ++i) {
						real rho = X[i] * X[i] + X[max_cnt + i] * X[max_cnt + i]
								+ X[2 * max_cnt + i] * X[2 * max_cnt + i];
						rho = std::sqrt(rho);
						real theta = std::acos(X[2 * max_cnt + i] / rho);
						real phi = std::atan2(X[1 * max_cnt + i], X[0 * max_cnt + i]);
						real x = std::cos(theta);                              // x = cos(theta)
						real y = std::sin(theta);                              // y = sin(theta)
						real fact = 1;                                   // Initialize 2 * m + 1
						real pn = 1;                        // Initialize Legendre polynomial Pn
						real rhom = 1.0 / rho;                          // Initialize rho^(-m-1)
#pragma novector
						for (int m = 0; m != P; ++m) {                     // Loop over m in Ynm
							real eim_r = std::cos(real(m) * phi);
							real eim_i = std::sin(real(m) * phi);
							real p = pn;                  //  Associated Legendre polynomial Pnm
							int npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
							int nmn = m * m;                          //  Index of Ynm for m < 0
							Ynm[npn * max_cnt + i] = rhom * p * _prefactor[npn] * eim_r; //  rho^(-m-1) * Ynm for m > 0
							if (npn != nmn) {
								Ynm[nmn * max_cnt + i] = rhom * p * _prefactor[npn] * eim_i; //  rho^(-m-1) * Ynm for m > 0
							}
							real p1 = p;                                              //  Pnm-1
							p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
							rhom /= rho;                                          //  rho^(-m-1)
							real rhon = rhom;                                     //  rho^(-n-1)
#pragma novector
							for (int n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
								int npm = n * n + n + m;             //   Index of Ynm for m > 0
								int nmm = n * n + n - m;             //   Index of Ynm for m < 0
								Ynm[npm * max_cnt + i] = rhon * p * _prefactor[npm] * eim_r; //   rho^n * Ynm for m > 0
								if (npm != nmm) {
									Ynm[nmm * max_cnt + i] = rhon * p * _prefactor[npm] * eim_i; //   rho^n * Ynm for m > 0
								}
								real p2 = p1;                                         //   Pnm-2
								p1 = p;                                               //   Pnm-1
								p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
								rhon /= rho;                                     //   rho^(-n-1)
							}                                         //  End loop over n in Ynm
							pn = -pn * fact * y;                                      //  Pn
							fact += 2;                                             //  2 * m + 1
						}                                              // End loop over m in Ynm
					}

					for (std::int64_t j = 0; j != P; ++j) {
						for (std::int64_t k = 0; k <= j; ++k) {
							const std::int64_t jkp = j * j + j + k;
							const std::int64_t jkm = j * j + j - k;
#pragma vector aligned
#pragma simd
							for (std::size_t i = 0; i != cnt; ++i) {
								L_r[i] = L_i[i] = real(0.0);
							}
							for (std::int64_t n = 0; n != P - j; ++n) {
								for (std::int64_t m = -n; m <= +n; ++m) {
									const std::int64_t nn = n * n + n;
									const std::int64_t nj = (n + j) * ((n + j) + 1);
									const std::int64_t jknm = jkp * P2 + n * n + n + m;
									const std::int64_t nmp = nn + std::abs(m);
									const std::int64_t nmm = nn - std::abs(m);
									const std::int64_t jnkmp = nj + std::abs(m - k);
									const std::int64_t jnkmm = nj - std::abs(m - k);
									real tmp_r, tmp_i;
									const real sgn = SGN(m-k);
									COMPLEX_MULT(tmp_r, tmp_i, M[nmp][j0 + Nx * (k0 + Nx * l0)],
											SGN(m) * M[nmm][j0 + Nx * (k0 + Nx * l0)], _Cnm_r[jknm], _Cnm_i[jknm]);
									const auto* Yp = Ynm + jnkmp * max_cnt;
									const auto* Ym = Ynm + jnkmm * max_cnt;
#pragma vector aligned
#pragma simd
									for (std::size_t i = 0; i != cnt; ++i) {
										COMPLEX_MULT_ADD(L_r[i], L_i[i], tmp_r, tmp_i, Yp[i], sgn * Ym[i]);
									}
								}
							}
#pragma vector aligned
#pragma simd
							for (std::size_t i = 0; i != cnt; ++i) {
								l_buffer[jkp * sz + indexes[i]] += L_r[i];
								l_buffer[jkm * sz + indexes[i]] += (k == 0) ? real(0.0) : L_i[i];
							}
						}
					}
				}
			}
		}
		delete[] buffer;
		delete[] indexes;
	}

	delete[] m_buffer;
	cnt = 0;
	for (std::size_t p = 0; p != P2; ++p) {
		for (std::size_t i = 0; i != sz; ++i) {
			L[p][i] = l_buffer[cnt++];
		}
	}
	delete[] l_buffer;
	return L;
}

template<std::int64_t P>
void exafmm_kernel<P>::M2L(std::valarray<std::valarray<real>>& CiL, const std::valarray<real> CjM,
		const std::valarray<std::valarray<real>>& d, std::size_t N) {
	std::valarray<std::valarray<real>> Ynm(std::valarray<real>(N), P * P);
#pragma vector aligned
#pragma simd
	for (std::size_t i = 0; i != N; ++i) {
		real rho = std::sqrt(d[0][i] * d[0][i] + d[1][i] * d[1][i] + d[2][i] * d[2][i]);
		real theta = std::acos(d[2][i] / rho);
		real phi = std::atan2(d[1][i], d[0][i]);
		real x = std::cos(theta);                              // x = cos(theta)
		real y = std::sin(theta);                              // y = sin(theta)
		real fact = 1;                                   // Initialize 2 * m + 1
		real pn = 1;                        // Initialize Legendre polynomial Pn
		real rhom = 1.0 / rho;                          // Initialize rho^(-m-1)
#pragma novector
		for (int m = 0; m != P; ++m) {                     // Loop over m in Ynm
			real eim_r = std::cos(real(m) * phi);
			real eim_i = std::sin(real(m) * phi);
			real p = pn;                  //  Associated Legendre polynomial Pnm
			int npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
			int nmn = m * m;                          //  Index of Ynm for m < 0
			Ynm[npn][i] = rhom * p * prefactor[npn] * eim_r; //  rho^(-m-1) * Ynm for m > 0
			if (npn != nmn) {
				Ynm[nmn][i] = rhom * p * prefactor[npn] * eim_i; //  rho^(-m-1) * Ynm for m > 0
			}
			real p1 = p;                                              //  Pnm-1
			p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
			rhom /= rho;                                          //  rho^(-m-1)
			real rhon = rhom;                                     //  rho^(-n-1)
#pragma novector
			for (int n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
				int npm = n * n + n + m;             //   Index of Ynm for m > 0
				int nmm = n * n + n - m;             //   Index of Ynm for m < 0
				Ynm[npm][i] = rhon * p * prefactor[npm] * eim_r; //   rho^n * Ynm for m > 0
				if (npm != nmm) {
					Ynm[nmm][i] = rhon * p * prefactor[npm] * eim_i; //   rho^n * Ynm for m > 0
				}
				real p2 = p1;                                         //   Pnm-2
				p1 = p;                                               //   Pnm-1
				p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
				rhon /= rho;                                     //   rho^(-n-1)
			}                                         //  End loop over n in Ynm
			pn = -pn * fact * y;                                      //  Pn
			fact += 2;                                             //  2 * m + 1
		}                                              // End loop over m in Ynm
	}

	std::valarray<real> L_r(N);
	std::valarray<real> L_i(N);
	for (std::int64_t j = 0; j != P; ++j) {
		for (std::int64_t k = 0; k <= j; ++k) {
			const std::int64_t jkp = j * j + j + k;
			const std::int64_t jkm = j * j + j - k;
#pragma vector aligned
#pragma simd
			for (std::size_t i = 0; i != N; ++i) {
				L_r[i] = L_i[i] = real(0.0);
			}
			for (std::int64_t n = 0; n != P - j; ++n) {
				for (std::int64_t m = -n; m <= +n; ++m) {
					const std::int64_t nn = n * n + n;
					const std::int64_t nj = (n + j) * ((n + j) + 1);
					const std::int64_t jknm = jkp * P * P + n * n + n + m;
					const std::int64_t nmp = nn + std::abs(m);
					const std::int64_t nmm = nn - std::abs(m);
					const std::int64_t jnkmp = nj + std::abs(m - k);
					const std::int64_t jnkmm = nj - std::abs(m - k);
					real tmp_r, tmp_i;
					const real sgn = SGN(m-k);
					COMPLEX_MULT(tmp_r, tmp_i, CjM[nmp], SGN(m) * CjM[nmm], Cnm_r[jknm], Cnm_i[jknm]);
					const auto& Yp = Ynm[jnkmp];
					const auto& Ym = Ynm[jnkmm];
#pragma vector aligned
#pragma simd
					for (std::size_t i = 0; i != N; ++i) {
						COMPLEX_MULT_ADD(L_r[i], L_i[i], tmp_r, tmp_i, Yp[i], sgn * Ym[i]);
					}
				}
			}
			auto& Cp = CiL[jkp];
			auto& Cm = CiL[jkm];
#pragma vector aligned
#pragma simd
			for (std::size_t i = 0; i != N; ++i) {
				Cp[i] = L_r[i];
				Cm[i] = (k == 0) ? L_r[i] : L_i[i];
			}
		}
	}

}

template<std::int64_t P>
void exafmm_kernel<P>::cart2sph(real& r, real& theta, real& phi, std::valarray<real> dist) {
	r = sqrt((dist * dist).sum()) * (1.0);      // r = sqrt(x^2 + y^2 + z^2)
	if (r < EPS) {                                             // If r == 0
		theta = 0;               //  theta can be anything so we set it to 0
	} else {                                                    // If r != 0
		theta = acos(dist[2] / r);                   //  theta = acos(z / r)
	}                                                   // End if for r == 0
	phi = atan2(dist[1], dist[0]);
}

template<std::int64_t P>
void exafmm_kernel<P>::M2M(std::valarray<std::valarray<real>>& CiM, const std::valarray<std::valarray<real>>& CjM,
		const std::valarray<real>& dist, const std::size_t N) {
	std::valarray<real> Ynm(P * P);
	std::valarray<real> M_r(N), M_i(N);
	real rho, theta, phi;
	cart2sph(rho, theta, phi, dist);
	evalMultipole(rho, theta, -phi, Ynm);
	for (std::int64_t j = 0; j != P; ++j) {
		for (std::int64_t k = 0; k <= j; ++k) {
			const std::int64_t jkp = j * j + j + k;
			const std::int64_t jkm = j * j + j - k;
#pragma vector aligned
#pragma simd
			for (std::size_t i = 0; i != N; ++i) {
				M_r[i] = M_i[i] = real(0.0);
			}
			for (std::int64_t n = 0; n <= j; ++n) {
				for (std::int64_t m = std::max(n - j + k, -n); m <= std::min(j - n + k, +n); ++m) {
					const std::int64_t nn = n * n + n;
					const std::int64_t nj = (j - n) * (j - n) + j - n;
					const std::int64_t jnkm = nj + k - m;
					const std::int64_t jnpkm = nj + std::abs(k - m);
					const std::int64_t jnmkm = nj - std::abs(k - m);
					const std::int64_t nmp = nn + std::abs(m);
					const std::int64_t nmm = nn - std::abs(m);
					const auto& Mj_r = CjM[jnpkm];
					const auto& Mj_i = CjM[jnmkm];
					const real tmp = Anm[nmp] * Anm[jnkm]
							/ Anm[jkp]* ODDEVEN((std::abs(k) - std::abs(m) - std::abs(k - m)) / 2 + n);
					const real sgn_km = SGN(k-m);
					const real Y_r = tmp * Ynm[nmp];
					const real Y_i = SGN(m) * tmp * Ynm[nmm];
#pragma vector aligned
#pragma simd
					for (std::size_t i = 0; i != N; ++i) {
						COMPLEX_MULT_ADD(M_r[i], M_i[i], Y_r, Y_i, Mj_r[i], sgn_km * Mj_i[i]);
					}
				}
			}
			auto& Mi_r = CiM[jkp];
			auto& Mi_i = CiM[jkm];
#pragma vector aligned
#pragma simd
			for (std::size_t i = 0; i != N; ++i) {
				Mi_r[i] = M_r[i];
				Mi_i[i] = (jkm == jkp) ? M_r[i] : M_i[i];
			}
		}
	}
}

template<std::int64_t P>
void exafmm_kernel<P>::L2L(std::valarray<std::valarray<real>>& CiL, const std::valarray<std::valarray<real>>& CjL,
		const std::valarray<real>& dist, const std::size_t N) {
	std::valarray<real> Ynm(P * P);
	real rho, theta, phi;
	std::valarray<real> L_r(N), L_i(N);
	cart2sph(rho, theta, phi, dist);
	evalMultipole(rho, theta, phi, Ynm);
	for (std::int64_t j = 0; j != P; ++j) {
		for (std::int64_t k = 0; k <= j; ++k) {
			std::int64_t jkp = j * j + j + k;
			std::int64_t jkm = j * j + j - k;
#pragma vector aligned
#pragma simd
			for (std::size_t i = 0; i != N; ++i) {
				L_r[i] = L_i[i] = 0.0;
			}
			for (std::int64_t n = j; n != P; ++n) {
				for (std::int64_t m = j - n + k; m <= n - j + k; ++m) {
					const std::int64_t nn = n * n + n;
					const std::int64_t nj = (n - j) * ((n - j) + 1);
					const std::int64_t npm = nn + std::abs(m);
					const std::int64_t nmm = nn - std::abs(m);
					const std::int64_t jnpkm = nj + std::abs(m - k);
					const std::int64_t jnmkm = nj - std::abs(m - k);
					const auto& Lj_r = CjL[npm];
					const auto& Lj_i = CjL[nmm];
					const real sgn = SGN(m);
					real tmp = std::pow(-1.0, (std::abs(m) - std::abs(k) - std::abs(m - k)) / 2) * Anm[jnpkm] * Anm[jkp]
							/ Anm[npm];
					const real Y_r = Ynm[jnpkm] * tmp;
					const real Y_i = SGN(m-k) * Ynm[jnmkm] * tmp;
#pragma vector aligned
#pragma simd
					for (std::size_t i = 0; i != N; ++i) {
						COMPLEX_MULT_ADD(L_r[i], L_i[i], Y_r, Y_i, Lj_r[i], sgn * Lj_i[i]);
					}
				}
			}
			auto& Li_r = CiL[jkp];
			auto& Li_i = CiL[jkm];
#pragma vector aligned
#pragma simd
			for (std::size_t i = 0; i != N; ++i) {
				Li_r[i] = L_r[i];
				Li_i[i] = (k == 0) ? L_r[i] : L_i[i];
			}
		}
	}
}

template<std::int64_t P>
void exafmm_kernel<P>::evalMultipole(real rho, real theta, real phi, std::valarray<real>& Ynm) {
	real x = std::cos(theta);                              // x = cos(theta)
	real y = std::sin(theta);                              // y = sin(theta)
	real fact = 1;                                   // Initialize 2 * m + 1
	real pn = 1;                        // Initialize Legendre polynomial Pn
	real rhom = 1;                                       // Initialize rho^m
	for (int m = 0; m != P; ++m) {                     // Loop over m in Ynm
		real eim_r = std::cos(real(m) * phi);
		real eim_i = std::sin(real(m) * phi);
		real p = pn;                  //  Associated Legendre polynomial Pnm
		int npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
		int nmn = m * m;                          //  Index of Ynm for m < 0
		Ynm[npn] = rhom * p * prefactor[npn] * eim_r; //  rho^m * Ynm for m > 0
		if (npn != nmn) {
			Ynm[nmn] = rhom * p * prefactor[npn] * eim_i; //  rho^m * Ynm for m > 0
		}
		real p1 = p;                                              //  Pnm-1
		p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
		rhom *= rho;                                              //  rho^m
		real rhon = rhom;                                         //  rho^n
		for (int n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
			int npm = n * n + n + m;             //   Index of Ynm for m > 0
			int nmm = n * n + n - m;             //   Index of Ynm for m < 0
			Ynm[npm] = rhon * p * prefactor[npm] * eim_r;     //   rho^n * Ynm
			if (npm != nmm) {
				Ynm[nmm] = rhon * p * prefactor[npm] * eim_i;     //   rho^n * Ynm
			}
			real p2 = p1;                                         //   Pnm-2
			p1 = p;                                               //   Pnm-1
			p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1); //   Pnm using recurrence relation
			rhon *= rho;                                   //   Update rho^n
		}                                         //  End loop over n in Ynm
		pn = -pn * fact * y;                                      //  Pn
		fact += 2;                                             //  2 * m + 1
	}                                              // End loop over m in Ynm
}

template<std::int64_t P>
exafmm_kernel<P>::exafmm_kernel() {
	const complex I(0., 1.);                               // Imaginary unit

	factorial[0] = 1;                                // Initialize factorial
	for (int n = 1; n != P; ++n) {                              // Loop to P
		factorial[n] = factorial[n - 1] * n;                        //  n!
	}                                                       // End loop to P

	for (int n = 0; n != P; ++n) {                     // Loop over n in Anm
		for (int m = -n; m <= n; ++m) {               //  Loop over m in Anm
			int nm = n * n + n + m;                        //   Index of Anm
			int nabsm = abs(m);                                     //   |m|
			real fnmm = 1.0;                        //   Initialize (n - m)!
			for (int i = 1; i <= n - m; ++i)
				fnmm *= i;                  //   (n - m)!
			real fnpm = 1.0;                        //   Initialize (n + m)!
			for (int i = 1; i <= n + m; ++i)
				fnpm *= i;                  //   (n + m)!
			real fnma = 1.0;                      //   Initialize (n - |m|)!
			for (int i = 1; i <= n - nabsm; ++i)
				fnma *= i;              //   (n - |m|)!
			real fnpa = 1.0;                      //   Initialize (n + |m|)!
			for (int i = 1; i <= n + nabsm; ++i)
				fnpa *= i;              //   (n + |m|)!
			prefactor[nm] = std::sqrt(fnma / fnpa); //   sqrt( (n - |m|)! / (n + |m|)! )
			Anm[nm] = ODDEVEN(n) / std::sqrt(fnmm * fnpm); //   (-1)^n / sqrt( (n + m)! / (n - m)! )
		}                                         //  End loop over m in Anm
	}                                              // End loop over n in Anm

	for (int j = 0, jk = 0, jknm = 0; j != P; ++j) { // Loop over j in Cjknm
		for (int k = -j; k <= j; ++k, ++jk) {       //  Loop over k in Cjknm
			for (int n = 0, nm = 0; n != P; ++n) { //   Loop over n in Cjknm
				for (int m = -n; m <= n; ++m, ++nm, ++jknm) { //    Loop over m in Cjknm
					if (j + n < P) {
						const int jnkm = (j + n) * (j + n) + j + n + m - k; //     Index C_{j+n}^{m-k}
						auto tmp = std::pow(I, real(abs(k - m) - abs(k) - abs(m))) //     Cjknm
						* real(ODDEVEN(j) * Anm[nm] * Anm[jk] / Anm[jnkm]);
						Cnm_r[jknm] = tmp.real();
						Cnm_i[jknm] = tmp.imag();

					}
				}
			}
		}
	}
}
