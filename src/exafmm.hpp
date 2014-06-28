/****
 * THIS FILE IS ADAPTED FROM EXAFMM (SEE COPYRIGHT NOTICE BELOW). THE INTERFACE HAS BEEN ALTERED
 * SLIGHTLY
 */

/*
 Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#include <valarray>
#ifndef kernel_h
#define kernel_h
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

#include <complex>

#define EPS2 0.0

using real = double;
using complex = std::complex<real>;

#define EPS (1.0e-12)

template<std::int64_t P>
class exafmm_kernel {
public:
	static void cart2sph(real& r, real& theta, real& phi, std::valarray<real> dist) {
		r = sqrt((dist * dist).sum()) * (1.0);                           // r = sqrt(x^2 + y^2 + z^2)
		if (r < EPS) {                                             // If r == 0
			theta = 0;                                                //  theta can be anything so we set it to 0
		} else {                                                    // If r != 0
			theta = acos(dist[2] / r);                                //  theta = acos(z / r)
		}                                                           // End if for r == 0
		phi = atan2(dist[1], dist[0]);
		/*	if (fabs(dist[0]) + fabs(dist[1]) < EPS) {                 // If |x| < eps & |y| < eps
		 phi = 0;                                                  //  phi can be anything so we set it to 0
		 } else if (fabs(dist[0]) < EPS) {                          // If |x| < eps
		 phi = dist[1] / fabs(dist[1]) * M_PI * 0.5;               //  phi = sign(y) * pi / 2
		 } else if (dist[0] > 0) {                                  // If x > 0
		 phi = atan(dist[1] / dist[0]);                            //  phi = atan(y / x)
		 } else {                                                    // If x < 0
		 phi = atan(dist[1] / dist[0]) + M_PI;                     //  phi = atan(y / x) + pi
		 }     */                                                      // End if for x,y cases
	}
	static void sph2cart(real r, real theta, real phi, std::valarray<real> spherical, std::valarray<real> &cartesian) {
		cartesian[0] = sin(theta) * cos(phi) * spherical[0]         // x component (not x itself)
		+ cos(theta) * cos(phi) / r * spherical[1] - sin(phi) / r / sin(theta) * spherical[2];
		cartesian[1] = sin(theta) * sin(phi) * spherical[0]         // y component (not y itself)
		+ cos(theta) * sin(phi) / r * spherical[1] + cos(phi) / r / sin(theta) * spherical[2];
		cartesian[2] = cos(theta) * spherical[0]                    // z component (not z itself)
		- sin(theta) / r * spherical[1];
	}

	static void M2M(std::valarray<complex>& CiM, const std::valarray<complex>& CjM, const std::valarray<real>& dist) {
		const complex I(0., 1.);
		complex Ynm[P * P], YnmTheta[P * P];

		real rho, alpha, beta;
		cart2sph(rho, alpha, beta, dist);
		evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
		for (int j = 0; j != P; ++j) {
			for (int k = 0; k <= j; ++k) {
				int jk = j * j + j + k;
				int jks = j * (j + 1) / 2 + k;
				complex M = 0;
				for (int n = 0; n <= j; ++n) {
					for (int m = -n; m <= std::min(k - 1, n); ++m) {
						if (j - n >= k - m) {
							int jnkm = (j - n) * (j - n) + j - n + k - m;
							int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
							int nm = n * n + n + m;
							M += CjM[jnkms] * std::pow(I, real(m - abs(m))) * Ynm[nm] * real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
						}
					}
					for (int m = k; m <= n; ++m) {
						if (j - n >= m - k) {
							int jnkm = (j - n) * (j - n) + j - n + k - m;
							int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
							int nm = n * n + n + m;
							M += std::conj(CjM[jnkms]) * Ynm[nm] * real(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
						}
					}
				}
				CiM[jks] += M;
			}
		}
	}

	static void M2L(std::valarray<complex> CiL, const std::valarray<complex>& CjM, const std::valarray<real>& dist) {
		complex Ynm[P * P], YnmTheta[P * P];
		real rho, alpha, beta;
		cart2sph(rho, alpha, beta, dist);
		evalLocal(rho, alpha, beta, Ynm, YnmTheta);
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
			}
		}
	}

	static void L2L(std::valarray<complex>& CiL, const std::valarray<complex>& CjL, const std::valarray<real>& dist) {
		const complex I(0., 1.);
		complex Ynm[P * P], YnmTheta[P * P];
		real rho, alpha, beta;
		cart2sph(rho, alpha, beta, dist);
		evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
		for (int j = 0; j != P; ++j) {
			for (int k = 0; k <= j; ++k) {
				int jk = j * j + j + k;
				int jks = j * (j + 1) / 2 + k;
				complex L = 0;
				for (int n = j; n != P; ++n) {
					for (int m = j + k - n; m < 0; ++m) {
						int jnkm = (n - j) * (n - j) + n - j + m - k;
						int nm = n * n + n - m;
						int nms = n * (n + 1) / 2 - m;
						L += std::conj(CjL[nms]) * Ynm[jnkm] * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
					}
					for (int m = 0; m <= n; ++m) {
						if (n - j >= abs(m - k)) {
							int jnkm = (n - j) * (n - j) + n - j + m - k;
							int nm = n * n + n + m;
							int nms = n * (n + 1) / 2 + m;
							L += CjL[nms] * std::pow(I, real(m - k - abs(m - k))) * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
						}
					}
				}
				CiL[jks] = L;
			}
		}
	}

	static void M2P(real& P0, std::valarray<real>& F0, const std::valarray<complex>& CjM, const std::valarray<real>& dist) {
		const complex I(0., 1.);                                       // Imaginary unit
		complex Ynm[P * P], YnmTheta[P * P];
		std::valarray<real> spherical(3, 0.0);
		std::valarray<real> cartesian(3, 0.0);
		real r, theta, phi;
		cart2sph(r, theta, phi, dist);
		evalLocal(r, theta, phi, Ynm, YnmTheta);

		for (int n = 0; n != P; ++n) {
			int nm = n * n + n;
			int nms = n * (n + 1) / 2;
			P0 += std::real(CjM[nms] * Ynm[nm]);
			spherical[0] -= std::real(CjM[nms] * Ynm[nm]) / r * (n + 1);
			spherical[1] += std::real(CjM[nms] * YnmTheta[nm]);
			for (int m = 1; m <= n; ++m) {
				nm = n * n + n + m;
				nms = n * (n + 1) / 2 + m;
				P0 += 2 * std::real(CjM[nms] * Ynm[nm]);
				spherical[0] -= 2 * std::real(CjM[nms] * Ynm[nm]) / r * (n + 1);
				spherical[1] += 2 * std::real(CjM[nms] * YnmTheta[nm]);
				spherical[2] += 2 * std::real(CjM[nms] * Ynm[nm] * I) * m;
			}
		}

		sph2cart(r, theta, phi, spherical, cartesian);
		F0 += cartesian;
	}
	static void L2P(real& P0, std::valarray<real>& F0, const std::valarray<complex>& CiL, const std::valarray<real>& dist) {
		const complex I(0., 1.);                                       // Imaginary unit
		complex Ynm[P * P], YnmTheta[P * P];
		std::valarray<real> spherical(3, 0.0);
		std::valarray<real> cartesian(3, 0.0);
		real r, theta, phi;
		cart2sph(r, theta, phi, dist);
		evalMultipole(r, theta, phi, Ynm, YnmTheta);
		for (int n = 0; n != P; ++n) {
			int nm = n * n + n;
			int nms = n * (n + 1) / 2;
			P0 += std::real(CiL[nms] * Ynm[nm]);
			spherical[0] += std::real(CiL[nms] * Ynm[nm]) / r * n;
			spherical[1] += std::real(CiL[nms] * YnmTheta[nm]);
			for (int m = 1; m <= n; ++m) {
				nm = n * n + n + m;
				nms = n * (n + 1) / 2 + m;
				P0 += 2 * std::real(CiL[nms] * Ynm[nm]);
				spherical[0] += 2 * std::real(CiL[nms] * Ynm[nm]) / r * n;
				spherical[1] += 2 * std::real(CiL[nms] * YnmTheta[nm]);
				spherical[2] += 2 * std::real(CiL[nms] * Ynm[nm] * I) * m;
			}
		}
		sph2cart(r, theta, phi, spherical, cartesian);
		F0 += cartesian;
	}

	static void P2P(real& P0, std::valarray<real>& F0, real m, std::valarray<real> dist) {         // Laplace P2P kernel on CPU
		real R2 = std::sqrt((dist * dist).sum());                              //   R^2
		real invR2 = 1.0 / R2;                                    //   1 / R^2
		if (R2 == 0) {											//   Exclude self interaction
			invR2 = 0;
		}
		real invR = m * std::sqrt(invR2);                   //   potential
		dist *= invR2 * invR;                                     //   force
		P0 += invR;                                               //   accumulate potential
		F0 += dist;                                               //   accumulate force
	}                                                             // End loop over target bodies

	static void evalMultipole(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) {
		const complex I(0., 1.);                                     // Imaginary unit
		real x = std::cos(alpha);                                   // x = cos(alpha)
		real y = std::sin(alpha);                                   // y = sin(alpha)
		real fact = 1;                                              // Initialize 2 * m + 1
		real pn = 1;                                                // Initialize Legendre polynomial Pn
		real rhom = 1;                                              // Initialize rho^m
		for (int m = 0; m != P; ++m) {                                 // Loop over m in Ynm
			complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
			real p = pn;                                              //  Associated Legendre polynomial Pnm
			int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
			int nmn = m * m;                                          //  Index of Ynm for m < 0
			Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^m * Ynm for m > 0
			Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
			real p1 = p;                                              //  Pnm-1
			p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
			YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;                                 // theta derivative of r^n * Ynm
			rhom *= rho;                                              //  rho^m
			real rhon = rhom;                                         //  rho^n
			for (int n = m + 1; n != P; ++n) {                             //  Loop over n in Ynm
				int npm = n * n + n + m;                                //   Index of Ynm for m > 0
				int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
				Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm
				Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
				real p2 = p1;                                           //   Pnm-2
				p1 = p;                                                 //   Pnm-1
				p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);                                                 //   Pnm using recurrence relation
				YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;                                      // theta derivative
				rhon *= rho;                                            //   Update rho^n
			}                                                         //  End loop over n in Ynm
			pn = -pn * fact * y;                                      //  Pn
			fact += 2;                                                //  2 * m + 1
		}                                                           // End loop over m in Ynm
	}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
	static void evalLocal(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) {
		const complex I(0., 1.);                                     // Imaginary unit
		real x = std::cos(alpha);                                   // x = cos(alpha)
		real y = std::sin(alpha);                                   // y = sin(alpha)
		real fact = 1;                                              // Initialize 2 * m + 1
		real pn = 1;                                                // Initialize Legendre polynomial Pn
		real rhom = 1.0 / rho;                                      // Initialize rho^(-m-1)
		for (int m = 0; m != P; ++m) {                                 // Loop over m in Ynm
			complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
			real p = pn;                                              //  Associated Legendre polynomial Pnm
			int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
			int nmn = m * m;                                          //  Index of Ynm for m < 0
			Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^(-m-1) * Ynm for m > 0
			Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
			real p1 = p;                                              //  Pnm-1
			p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
			YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;                                 // theta derivative of r^n * Ynm
			rhom /= rho;                                              //  rho^(-m-1)
			real rhon = rhom;                                         //  rho^(-n-1)
			for (int n = m + 1; n != P; ++n) {                             //  Loop over n in Ynm
				int npm = n * n + n + m;                                //   Index of Ynm for m > 0
				int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
				Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm for m > 0
				Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
				real p2 = p1;                                           //   Pnm-2
				p1 = p;                                                 //   Pnm-1
				p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);                                                 //   Pnm using recurrence relation
				YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;                                      // theta derivative
				rhon /= rho;                                            //   rho^(-n-1)
			}                                                         //  End loop over n in Ynm
			pn = -pn * fact * y;                                      //  Pn
			fact += 2;                                                //  2 * m + 1
		}                                                           // End loop over m in Ynm
	}

private:

	static std::array<real, P> factorial;
	static std::array<real, P * P> prefactor;
	static std::array<real, P * P> Anm;
	static std::array<complex, P * P * P * P> Cnm;

public:
	exafmm_kernel() {
		const complex I(0., 1.);                                     // Imaginary unit

		factorial[0] = 1;                                           // Initialize factorial
		for (int n = 1; n != P; ++n) {                                 // Loop to P
			factorial[n] = factorial[n - 1] * n;                        //  n!
		}                                                           // End loop to P

		for (int n = 0; n != P; ++n) {                                 // Loop over n in Anm
			for (int m = -n; m <= n; ++m) {                              //  Loop over m in Anm
				int nm = n * n + n + m;                                       //   Index of Anm
				int nabsm = abs(m);                                     //   |m|
				real fnmm = 1.0;                                        //   Initialize (n - m)!
				for (int i = 1; i <= n - m; ++i)
					fnmm *= i;                  //   (n - m)!
				real fnpm = 1.0;                                        //   Initialize (n + m)!
				for (int i = 1; i <= n + m; ++i)
					fnpm *= i;                  //   (n + m)!
				real fnma = 1.0;                                        //   Initialize (n - |m|)!
				for (int i = 1; i <= n - nabsm; ++i)
					fnma *= i;              //   (n - |m|)!
				real fnpa = 1.0;                                        //   Initialize (n + |m|)!
				for (int i = 1; i <= n + nabsm; ++i)
					fnpa *= i;              //   (n + |m|)!
				prefactor[nm] = std::sqrt(fnma / fnpa);                   //   sqrt( (n - |m|)! / (n + |m|)! )
				Anm[nm] = ODDEVEN(n) / std::sqrt(fnmm * fnpm);              //   (-1)^n / sqrt( (n + m)! / (n - m)! )
			}                                                         //  End loop over m in Anm
		}                                                           // End loop over n in Anm

		for (int j = 0, jk = 0, jknm = 0; j != P; ++j) {                   // Loop over j in Cjknm
			for (int k = -j; k <= j; ++k, ++jk) {                         //  Loop over k in Cjknm
				for (int n = 0, nm = 0; n != P; ++n) {                       //   Loop over n in Cjknm
					for (int m = -n; m <= n; ++m, ++nm, ++jknm) {            //    Loop over m in Cjknm
						if (j + n < P) {
							const int jnkm = (j + n) * (j + n) + j + n + m - k;               //     Index C_{j+n}^{m-k}
							Cnm[jknm] = std::pow(I, real(abs(k - m) - abs(k) - abs(m)))               //     Cjknm
							* real(ODDEVEN(j) * Anm[nm] * Anm[jk] / Anm[jnkm]);
						}                                                     //    End loop over m in Cjknm
					}
				}             //   End loop over n in Cjknm
			}                                                         //  End loop over in k in Cjknm
		}                                                           // End loop over in j in Cjknm
	}
};

template<std::int64_t P>
std::array<real, P> exafmm_kernel<P>::factorial;

template<std::int64_t P>
std::array<real, P * P> exafmm_kernel<P>::prefactor;

template<std::int64_t P>
std::array<real, P * P> exafmm_kernel<P>::Anm;

template<std::int64_t P>
std::array<complex, P * P * P * P> exafmm_kernel<P>::Cnm;

#endif
