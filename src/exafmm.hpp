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
#include <array>

#define EPS2 0.0

using real = double;
using complex = std::complex<real>;

#define EPS (1.0e-12)

template<std::int64_t P>
class exafmm_kernel {
public:
	static void cart2sph(real& r, real& theta, real& phi, std::valarray<real> dist);
	static void M2M(std::valarray<real>& CiM, const std::valarray<real>& CjM, const std::valarray<real>& dist);
	static void M2L(std::valarray<real>& CiL, const std::valarray<real> CjM, const std::valarray<real>& dist);
	static void L2L(std::valarray<real>& CiL, const std::valarray<real>& CjL, const std::valarray<real>& dist);
	static void evalMultipole(real rho, real theta, real phi, std::valarray<real>& Ynm);
	static void evalLocal(real rho, real theta, real phi, std::valarray<real>& Ynm);
private:

	static std::array<real, P> factorial;
	static std::array<real, P * P> prefactor;
	static std::array<real, P * P> Anm;
	static std::array<complex, P * P * P * P> Cnm;

public:
	exafmm_kernel();
};

template<std::int64_t P>
std::array<real, P> exafmm_kernel<P>::factorial;

template<std::int64_t P>
std::array<real, P * P> exafmm_kernel<P>::prefactor;

template<std::int64_t P>
std::array<real, P * P> exafmm_kernel<P>::Anm;

template<std::int64_t P>
std::array<complex, P * P * P * P> exafmm_kernel<P>::Cnm;

#ifndef EXAFMM_CPP
extern template class exafmm_kernel<1>;
extern template class exafmm_kernel<2>;
extern template class exafmm_kernel<3>;
extern template class exafmm_kernel<4>;
extern template class exafmm_kernel<5>;
extern template class exafmm_kernel<6>;
extern template class exafmm_kernel<7>;
extern template class exafmm_kernel<8>;
extern template class exafmm_kernel<9>;
extern template class exafmm_kernel<10>;
#endif

#endif
