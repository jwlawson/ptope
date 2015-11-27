/*
 * duplicate_column_check.h
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
#pragma once
#ifndef PTOPE_DUPLICATE_COL_CHECK_H_
#define PTOPE_DUPLICATE_COL_CHECK_H_

#include "polytope_candidate.h"

namespace ptope {
/**
 * Check whether the given polytope has duplicate vectors. This is assumed to
 * be run at every step, so only the last added vector is checked to see if it
 * has already been added.
 *
 * Returns true if the polytope has duplicate vectors.
 */
class DuplicateColumnCheck {
static constexpr double error = 1e-14;
public:
	bool operator()(const PolytopeCandidate & p) {
		return operator()(p.gram());
	}
	bool operator()(const arma::mat & gram) {
		const arma::uword last_col_ind = gram.n_cols - 1;
		const arma::vec & last_col = gram.unsafe_col(last_col_ind);
		arma::vec diff(last_col.size());
		bool no_duplicate = true;
		for(arma::uword i = 0; no_duplicate && i < last_col_ind; ++i) {
			const arma::vec & col = gram.unsafe_col(i);
			diff = last_col - col;
			bool equal = true;
			for(const double & el : diff) {
				if(std::abs(el) > error) {
					equal = false;
					break;
				}
			}
			if(equal) {
				no_duplicate = false;
			}
		}
		return !no_duplicate;
	}
};
}
#endif

