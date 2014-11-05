/*
 * exafmm.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: dmarce1
 */

#define EXAFMM_CPP
#include "exafmm.hpp"

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
void exafmm_kernel<P>::M2M(std::valarray<real>& CiM, const std::valarray<real>& CjM, const std::valarray<real>& dist) {
	const complex I(0., 1.);
	std::valarray<real> Ynm(P * P);

	real rho, theta, phi;
	cart2sph(rho, theta, phi, dist);
	evalMultipole(rho, theta, -phi, Ynm);
	for (int j = 0; j != P; ++j) {
		for (int k = 0; k <= j; ++k) {
			int jk = j * j + j + k;
			int jks = j * (j + 1) / 2 + k;
			complex M = 0;
			for (int n = 0; n <= j; ++n) {
				for (int m = -n; m <= std::min(k - 1, n); ++m) {
					if (j - n >= k - m) {
						int jnkm = (j - n) * (j - n) + j - n + k - m;
						int nm = n * n + n + m;
						int jnpkm = (j - n) * ((j - n) + 1) + std::abs(k - m);
						int jnmkm = (j - n) * ((j - n) + 1) - std::abs(k - m);
						complex jM(CjM[jnpkm], +CjM[jnmkm]);
						auto Y = complex(Ynm[n * (n + 1) + std::abs(m)], m == 0 ? 0.0 : Ynm[n * (n + 1) - std::abs(m)]);
						if (m < 0) {
							Y = std::conj(Y);
						}
						M += jM * std::pow(I, real(m - abs(m))) * Y * real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
					}
				}
				for (int m = k; m <= n; ++m) {
					if (j - n >= m - k) {
						int jnkm = (j - n) * (j - n) + j - n + k - m;
						int nm = n * n + n + m;
						int jnpkm = (j - n) * ((j - n) + 1) + std::abs(k - m);
						int jnmkm = (j - n) * ((j - n) + 1) - std::abs(k - m);
						complex jM(CjM[jnpkm], k == m ? real(0.0) : -CjM[jnmkm]);
						auto Y = complex(Ynm[n * (n + 1) + std::abs(m)], m == 0 ? 0.0 : Ynm[n * (n + 1) - std::abs(m)]);
						if (m < 0) {
							Y = std::conj(Y);
						}
						M += jM * Y * real(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
					}
				}
			}
			CiM[j * j + j + k] += M.real();
			if (k != 0) {
				CiM[j * j + j - k] += M.imag();
			}
		}
	}
}

template<std::int64_t P>
void exafmm_kernel<P>::M2L(std::valarray<real>& CiL, const std::valarray<real> CjM, const std::valarray<real>& dist) {
	std::valarray<real> Ynm(P * P);
	real rho, theta, phi;
	cart2sph(rho, theta, phi, dist);
	evalLocal(rho, theta, phi, Ynm);
	for (int j = 0; j != P; ++j) {
		for (int k = 0; k <= j; ++k) {
			int jk = j * j + j + k;
			complex L = 0;
			for (int n = 0; n != P - j; ++n) {
				for (int m = -n; m < 0; ++m) {
					int nm = n * n + n + m;
					int nms = n * (n + 1) / 2 - m;
					int jknm = jk * P * P + nm;
					int jnkm = (j + n) * (j + n) + j + n + m - k;
					complex jM(CjM[n * n + n - m], -CjM[n * n + n + m]);
					auto Y = complex(Ynm[(n + j) * ((n + j) + 1) + std::abs(m - k)],
							m - k == 0 ? 0.0 : Ynm[(n + j) * ((n + j) + 1) - std::abs(m - k)]);
					if (m - k < 0) {
						Y = std::conj(Y);
					}
					L += jM * complex(Cnm_r[jknm], Cnm_i[jknm]) * Y;
				}
				for (int m = 0; m <= n; ++m) {
					int nm = n * n + n + m;
					int nms = n * (n + 1) / 2 + m;
					int jknm = jk * P * P + nm;
					int jnkm = (j + n) * (j + n) + j + n + m - k;
					complex jM(CjM[n * n + n + m], m == 0 ? real(0.0) : +CjM[n * n + n - m]);
					auto Y = complex(Ynm[(n + j) * ((n + j) + 1) + std::abs(m - k)],
							m - k == 0 ? 0.0 : Ynm[(n + j) * ((n + j) + 1) - std::abs(m - k)]);
					if (m - k < 0) {
						Y = std::conj(Y);
					}
					L += jM * complex(Cnm_r[jknm], Cnm_i[jknm]) * Y;
				}
			}
			CiL[j * j + j + k] = L.real();
			if (k != 0) {
				CiL[j * j + j - k] = L.imag();
			}
		}
	}
}

template<std::int64_t P>
void exafmm_kernel<P>::L2L(std::valarray<real>& CiL, const std::valarray<real>& CjL, const std::valarray<real>& dist) {
	const complex I(0., 1.);
	std::valarray<real> Ynm(P * P);
	real rho, theta, phi;
	cart2sph(rho, theta, phi, dist);
	evalMultipole(rho, theta, phi, Ynm);
	for (int j = 0; j != P; ++j) {
		for (int k = 0; k <= j; ++k) {
			int jk = j * j + j + k;
			complex L = 0;
			for (int n = j; n != P; ++n) {
				for (int m = j + k - n; m < 0; ++m) {
					complex jL(CjL[n * n + n - m], CjL[n * n + n + m]);
					int jnkm = (n - j) * (n - j) + n - j + m - k;
					int nm = n * n + n - m;
					auto Y = complex(Ynm[(n - j) * ((n - j) + 1) + std::abs(m - k)],
							m - k == 0 ? 0.0 : Ynm[(n - j) * ((n - j) + 1) - std::abs(m - k)]);
					if (m - k < 0) {
						Y = std::conj(Y);
					}
					L += std::conj(jL) * Y * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
				}
				for (int m = 0; m <= n; ++m) {
					if (n - j >= abs(m - k)) {
						complex jL(CjL[n * n + n + m], m == 0 ? 0.0 : CjL[n * n + n - m]);
						int jnkm = (n - j) * (n - j) + n - j + m - k;
						int nm = n * n + n + m;
						auto Y = complex(Ynm[(n - j) * ((n - j) + 1) + std::abs(m - k)],
								m - k == 0 ? 0.0 : Ynm[(n - j) * ((n - j) + 1) - std::abs(m - k)]);
						if (m - k < 0) {
							Y = std::conj(Y);
						}
						L += jL * std::pow(I, real(m - k - abs(m - k))) * Y * Anm[jnkm] * Anm[jk] / Anm[nm];
					}
				}
			}
			CiL[j * j + j + k] = L.real();
			if (k != 0) {
				CiL[j * j + j - k] = L.imag();
			}
		}
	}
}

template<std::int64_t P>
void exafmm_kernel<P>::evalMultipole(real rho, real theta, real phi, std::valarray<real>& Ynm) {
	const complex I(0., 1.);                               // Imaginary unit
	real x = std::cos(theta);                              // x = cos(theta)
	real y = std::sin(theta);                              // y = sin(theta)
	real fact = 1;                                   // Initialize 2 * m + 1
	real pn = 1;                        // Initialize Legendre polynomial Pn
	real rhom = 1;                                       // Initialize rho^m
	for (int m = 0; m != P; ++m) {                     // Loop over m in Ynm
		complex eim = std::exp(I * real(m * phi));      //  exp(i * m * phi)
		real p = pn;                  //  Associated Legendre polynomial Pnm
		int npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
		int nmn = m * m;                          //  Index of Ynm for m < 0
		Ynm[npn] = rhom * p * prefactor[npn] * eim.real(); //  rho^m * Ynm for m > 0
		if (npn != nmn) {
			Ynm[nmn] = rhom * p * prefactor[npn] * eim.imag(); //  rho^m * Ynm for m > 0
		}
		real p1 = p;                                              //  Pnm-1
		p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
		rhom *= rho;                                              //  rho^m
		real rhon = rhom;                                         //  rho^n
		for (int n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
			int npm = n * n + n + m;             //   Index of Ynm for m > 0
			int nmm = n * n + n - m;             //   Index of Ynm for m < 0
			Ynm[npm] = rhon * p * prefactor[npm] * eim.real();     //   rho^n * Ynm
			if (npm != nmm) {
				Ynm[nmm] = rhon * p * prefactor[npm] * eim.imag();     //   rho^n * Ynm
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

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
template<std::int64_t P>
void exafmm_kernel<P>::evalLocal(real rho, real theta, real phi, std::valarray<real>& Ynm) {
	const complex I(0., 1.);                               // Imaginary unit
	real x = std::cos(theta);                              // x = cos(theta)
	real y = std::sin(theta);                              // y = sin(theta)
	real fact = 1;                                   // Initialize 2 * m + 1
	real pn = 1;                        // Initialize Legendre polynomial Pn
	real rhom = 1.0 / rho;                          // Initialize rho^(-m-1)
	for (int m = 0; m != P; ++m) {                     // Loop over m in Ynm
		complex eim = std::exp(I * real(m * phi));      //  exp(i * m * phi)
		real p = pn;                  //  Associated Legendre polynomial Pnm
		int npn = m * m + 2 * m;                  //  Index of Ynm for m > 0
		int nmn = m * m;                          //  Index of Ynm for m < 0
		Ynm[npn] = rhom * p * prefactor[npn] * eim.real(); //  rho^(-m-1) * Ynm for m > 0
		if (npn != nmn) {
			Ynm[nmn] = rhom * p * prefactor[npn] * eim.imag(); //  rho^(-m-1) * Ynm for m > 0
		}
		real p1 = p;                                              //  Pnm-1
		p = x * (2 * m + 1) * p1;          //  Pnm using recurrence relation
		rhom /= rho;                                          //  rho^(-m-1)
		real rhon = rhom;                                     //  rho^(-n-1)
		for (int n = m + 1; n != P; ++n) {            //  Loop over n in Ynm
			int npm = n * n + n + m;             //   Index of Ynm for m > 0
			int nmm = n * n + n - m;             //   Index of Ynm for m < 0
			Ynm[npm] = rhon * p * prefactor[npm] * eim.real(); //   rho^n * Ynm for m > 0
			if (npm != nmm) {
				Ynm[nmm] = rhon * p * prefactor[npm] * eim.imag(); //   rho^n * Ynm for m > 0
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
