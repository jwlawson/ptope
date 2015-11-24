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
bool ParabolicCheck::operator()(const arma::mat & m) {
	if(std::abs(arma::det(m)) > error) {
		return false;
	}
	bool got_vals = arma::eig_sym(evalues, m);
	if(!got_vals) {
		std::cerr << "Failed to find eigenvalues: " << std::endl << m << std::endl;
	}
	bool is_non_negative = true;
	for(arma::uword i = 0; is_non_negative && i < m.n_rows; ++i) {
		if (evalues(i) < 0) {
			is_non_negative = false;
		}
	}
	return is_non_negative;
}
bool ParabolicCheck::operator()(const PolytopeCandidate & p) {
	return operator()(p.gram());
}
}
