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
bool ParabolicCheck::operator()(const arma::mat & m, const arma::uword & dim) {
	arma::uvec ind(dim);
	return check_submatrices(ind, 0, m);
}
bool ParabolicCheck::operator()(const PolytopeCandidate & p) {
	return operator()(p.gram(), p.real_dimension());
}
/*
 * Currently uses stupid eigendecomposition way of determining whether
 * parabolic. Should instead use something like Cholesky or LDL.
 */
bool
ParabolicCheck::parabolic(const arma::mat & m) {
	bool result = true;
	if(std::abs(arma::det(m)) > error) {
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
		result = parabolic(m.submat(indices, indices));
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
}
