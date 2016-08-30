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

#include "angles.h"

namespace ptope {
namespace {
constexpr double error = 1e-10;
}
AngleCheck::AngleCheck()
: _values(Angles::get().inner_products()) {}
bool
AngleCheck::operator()(const PolytopeCandidate & p) {
	return operator()(p.gram());
}
bool
AngleCheck::operator()(const arma::mat & m) {
	bool result = true;
	const arma::uword last_col_ind = m.n_cols - 1;
	const arma::vec & last_col = m.unsafe_col(last_col_ind);
	for(arma::uword i = 0; result && i < last_col_ind; ++i) {
		const double & val = last_col(i);
		result = operator()(val);
	}
	return result;
}
bool
AngleCheck::operator()(const double & val) {
	bool is_valid = 
		val < error + 0.0  // No entries apart from diagonals should be greater than 0
		&& (  val < error + -1.0   // No angle
			|| std::binary_search(_values.begin(), _values.end(), val, _dless) );
	return is_valid;
}
}
