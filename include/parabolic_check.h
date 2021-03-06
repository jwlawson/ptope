/*
 * parabolic_check.h
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
/**
 * Functor to check whether a given connected symmetric matrix contains a
 * parabolic submatrix including the last row/column.
 *
 * A connected matrix is parabolic if it is positive semi-definite but not
 * positive definite. Equivalently all eigenvalues are non-negative, with at
 * least one being 0.
 */
#pragma once
#ifndef PTOPE_PARABOLIC_CHECK_H_
#define PTOPE_PARABOLIC_CHECK_H_

#include "polytope_candidate.h"

namespace ptope {
class ParabolicCheck {
private:
	/**
	 * Class to compute determinants (up to sign) which keeps cetain chunks of
	 * memory allocated to prevent allocations each time a det is computed.
	 */
	class UnsignedDet {
	public:
		/**
		 * Compute the determinant of a given matrix. This determinant will give a
		 * correct absolute value, however the sign might be incorrect.
		 */
		double
		operator()(const arma::mat & m);
	private:
		arma::podarray<arma::blas_int> _ipiv;
		arma::mat _copy;
	};
public:
	ParabolicCheck() {}
	/**
	 * Check whether the given symmetric matrix contains a parabolic submatrix
	 * of the specified dimension including the last column/row.
	 *
	 * Returns true if the matrix does contain a parabolic submatrix.
	 */
	bool
	operator()(const arma::mat & m, const arma::uword & dim);
	bool
	operator()(const PolytopeCandidate & p);
private:
	arma::vec _evalues;
	UnsignedDet _det;
	arma::mat _submat_cache;
	/** Checks whether given matrix is parabolic. */
	bool
	parabolic(const arma::mat & m);
	/** Recursively construct submatrices. */
	bool
	check_submatrices(arma::uvec & indices, arma::uword index,
			const arma::mat & m);
};
inline
bool
ParabolicCheck::operator()(const arma::mat & m, const arma::uword & dim) {
	arma::uvec ind(dim);
	return check_submatrices(ind, 0, m);
}
inline
bool
ParabolicCheck::operator()(const PolytopeCandidate & p) {
	return operator()(p.gram(), p.real_dimension());
}
}
#endif

