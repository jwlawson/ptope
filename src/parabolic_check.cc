/*
 * parabolic_check.cc
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
#include "parabolic_check.h"

namespace ptope {
namespace {
constexpr double error = 1e-10;
}
/*
 * Currently uses stupid eigendecomposition way of determining whether
 * parabolic. Should instead use something like Cholesky or LDL.
 */
bool
ParabolicCheck::parabolic(const arma::mat & m) {
	bool result = true;
	if(std::abs(_det(m)) > error) {
		result = false;
	} else {
		arma::eig_sym(_evalues, m);
		for(arma::uword i = 0; result && i < m.n_rows; ++i) {
			if (_evalues(i) < -error) {
				result = false;
			}
		}
	}
	return result;
}
bool
ParabolicCheck::check_submatrices(arma::uvec & indices, arma::uword index,
		const arma::mat & m) {
	bool result = false;
	if(index == indices.size() - 1) {
		/* Have d-1 submatrix, so just add last column */
		indices(indices.size() - 1) = m.n_cols - 1;
		_submat_cache = m.submat(indices, indices);
		result = parabolic(_submat_cache);
	} else {
		/* Keep adding to submatrix */
		arma::uword min = (index == 0 ? 0 : indices(index - 1) + 1);
		for(arma::uword k = min, max = m.n_cols - 1; !result && k < max; ++k) {
			indices(index) = k;
			result = check_submatrices(indices, index + 1, m);
		}
	}
	return result;
}
double
ParabolicCheck::UnsignedDet::operator()(const arma::mat & m) {
	double val;
	if(m.n_rows <= 4) {
		val = arma::auxlib::det_tinymat(m, m.n_rows);
	} else {
		_copy = m;
		_ipiv.set_min_size(_copy.n_rows);
			
		arma::blas_int info(0);
		arma::blas_int n_rows(_copy.n_rows);
		arma::blas_int n_cols(_copy.n_cols);
			
		/* getrf computes LU decomposition. ipiv returns info on any row
		 * permutations (each of which would switch sign of det) */
		arma::lapack::getrf(&n_rows, &n_cols, _copy.memptr(), &n_rows,
				_ipiv.memptr(), &info);
		if(info > 0) {
			val = 0;
		} else {
			val = _copy.at(0,0);
			for(arma::uword i=1; i < _copy.n_rows; ++i){
				val *= _copy.at(i,i);
			}
		}
	}
	return val;
}
}
