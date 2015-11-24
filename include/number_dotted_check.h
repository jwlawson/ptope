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
#ifndef PTOPE_NUMBER_DOTED_CHECK_H_
#define PTOPE_NUMBER_DOTED_CHECK_H_

#include "polytope_candidate.h"

namespace ptope {
template <int N>
class NumberDottedCheck {
	public:
		bool operator()(const PolytopeCandidate & p) {
			const arma::mat & gram = p.gram();
			int count;
			for(const double & val : gram) {
				if(val > 1.0) {
					++count;
				}
			}
			return count == N;
		}
};
}
#endif

