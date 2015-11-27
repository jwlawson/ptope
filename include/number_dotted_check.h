/*
 * number_dotted_check.h
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
#ifndef PTOPE_NUMBER_DOTTED_CHECK_H_
#define PTOPE_NUMBER_DOTTED_CHECK_H_

#include "polytope_candidate.h"

namespace ptope {
/**
 * Checks the number of dotted edges in the last column of the gram matrix of
 * the polytope.
 */
template <int N>
class NumberDottedCheck {
	static constexpr double error = 1e-14;
	public:
		bool operator()(const PolytopeCandidate & p) {
			const arma::mat & gram = p.gram();
			const arma::uword last_col = gram.n_cols - 1;
			const arma::vec & col = gram.unsafe_col(last_col);
			int count;
			for(const double & val : col) {
				if(val - 1.0 > error) {
					++count;
				}
			}
			return count == N;
		}
};
}
#endif

