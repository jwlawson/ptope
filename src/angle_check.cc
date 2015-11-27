/*
 * angle_check.cc
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
#include "angle_check.h"

namespace ptope {
namespace {
constexpr double error = 1e-10;
std::vector<double> angles_to_prods(const std::vector<uint> & angles) {
	std::vector<double> result(angles.size() + 1);
	for(std::size_t i = 0, max = angles.size(); i < max; ++i) {
		if(angles[i] == 2) {
			result[i] = 0;
		} else {
			result[i] = -std::cos(arma::datum::pi/angles[i]);
		}
	}
	/* Always want an entry 1 to be valid. */
	result[angles.size()] = 1;
	return result;
}
}
const std::vector<uint>
AngleCheck::__default_mults = { 2, 3, 4, 5, 8 };
const std::vector<double>
AngleCheck::__default_inner = angles_to_prods(__default_mults);
AngleCheck::AngleCheck()
	: _values(__default_inner) {}
AngleCheck::AngleCheck(const std::vector<uint> & angles)
	: _values(angles_to_prods(angles)) {}
AngleCheck::AngleCheck(std::vector<uint> && angles)
	: _values(angles_to_prods(angles)) {}
bool
AngleCheck::operator()(const PolytopeCandidate & p) {
	return operator()(p.gram());
}
bool
AngleCheck::operator()(const arma::mat & m) {
	const arma::uword last_col_ind = m.n_cols - 1;
	const arma::vec & last_col = m.unsafe_col(last_col_ind);
	for(const double & val : last_col) {
		/* Uses lambda to check double values up to error */
		if( val < 1.0 && std::find_if(_values.begin(), _values.end(),
					[val](const double & b){return std::abs(b - val) < error; } )
				== _values.end()) {
			/* Value not found */
			return false;
		}
	}
	return true;
}
}
