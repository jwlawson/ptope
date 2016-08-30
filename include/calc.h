/*
 * calc.h
 * Copyright 2015 John Lawson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*
 * Provides calculations and simple functions which varioous classes need to
 * evaluate. For example min_cos_angle and mik_inner_product.
 */
#pragma once
#ifndef _PTOPE_CALC_H_
#define _PTOPE_CALC_H_

#include <armadillo>

namespace ptope {
namespace calc {
double
inline
min_cos_angle(unsigned int mult) {
	return -std::cos(arma::datum::pi / mult);
}
double
inline
eucl_inner_prod(std::size_t n_elems, double const * const a,
		double const * const b) {
	double sq1(0);
	double sq2(0);
	arma::uword i(0);
	for(arma::uword j = 1; j < n_elems; i+=2, j+=2) {
		sq1 += a[i] * b[i];
		sq2 += a[j] * b[j];
	}
	if(i < n_elems) {
		sq1 += a[i] * b[i];
	}
	return sq1 + sq2;
}
double
inline
eucl_inner_prod(const arma::vec & a, const arma::vec & b) {
	return std::inner_product(a.begin(), a.end(), b.begin(),
			static_cast<double>(0));
}
double
inline
eucl_sq_norm(const arma::vec & a) {
	return std::inner_product(a.begin(), a.end(), a.begin(),
			static_cast<double>(0));
}
/** Compute Minkowski/Lorentz inner product of two vectors. */
double
inline
mink_inner_prod(std::size_t n_elems, double const * const a,
		double const * const b) {
	const std::size_t max = n_elems - 1;
	double result = eucl_inner_prod(max, a, b);
	double m = a[max] * b[max];
	result -= m;
	return result;
}
double
inline
mink_inner_prod(const arma::vec & a, const arma::vec & b) {
	const std::size_t max = a.size() - 1;
	double result;
	switch(max){
		case 2:
			result = a[0]*b[0] + a[1]*b[1];
			break;
		case 3:
			result = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
			break;
		case 4:
			result = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
			break;
		case 5:
			result = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3] + a[4]*b[4];
			break;
		case 6:
			result = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3] + a[4]*b[4]
				+ a[5]*b[5];
			break;
		case 7:
			result = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3] + a[4]*b[4]
				+ a[5]*b[5] + a[6]*b[6];
			break;
		default:
			result = std::inner_product(a.begin(), a.end() - 1, b.begin(),
			static_cast<double>(0));
	}
	double m = a[max] * b[max];
	result -= m;
	return result;
}
double
inline
mink_sq_norm(const arma::vec & a) {
	const std::size_t max = a.size() - 1;
	double eucl = std::inner_product(a.begin(), a.end() - 1, a.begin(),
			static_cast<double>(0));
	double m = a[max] * a[max];
	return eucl - m;
}

}
}
#endif
