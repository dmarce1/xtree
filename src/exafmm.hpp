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
#include <hpx/runtime/actions/plain_action.hpp>
#include <valarray>
#include <boost/serialization/valarray.hpp>
#include <boost/serialization/complex.hpp>
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
	static std::valarray<std::valarray<complex>> M2L_interior(std::valarray<std::valarray<complex>> M, real dx, std::int64_t N, bool, bool );

	static void cart2sph(real& r, real& theta, real& phi, std::valarray<real> dist);
	static void M2M(std::valarray<complex>& CiM, const std::valarray<complex>& CjM, const std::valarray<real>& dist);
	static void M2M_vec(std::valarray<std::valarray<complex>>& CiM, const std::valarray<std::valarray<complex>>& CjM,
			real x0, real y0, real z0, std::size_t vlen);
	static void M2L(std::valarray<complex>& CiL, const std::valarray<complex> CjM, const std::valarray<real>& dist);
	static std::valarray<std::valarray<complex>> M2L_vec(const std::valarray<complex> CjM,
			const std::valarray<real>&, const std::valarray<real>&, const std::valarray<real>&, std::size_t vlens);
	static void L2L(std::valarray<complex>& CiL, const std::valarray<complex>& CjL, const std::valarray<real>& dist);
	static void evalMultipole(real rho, real theta, real phi, std::valarray<complex>& Ynm);
	static void evalLocal(real rho, real theta, real phi, std::valarray<complex>& Ynm);
private:

	static std::valarray<real> factorial;
	static std::valarray<real> prefactor;
	static std::valarray<real> Anm;
	static std::valarray<complex> Cnm;
	static std::valarray<real> Cnm_real;
	static std::valarray<real> Cnm_imag;

public:
	exafmm_kernel();
}
;

HPX_PLAIN_ACTION(exafmm_kernel<5>::M2L_interior, M2L_interior_action);

#endif

