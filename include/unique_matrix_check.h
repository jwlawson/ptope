/*
 * unique_matrix_check.h
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
#ifndef PTOPE_UNIQUE_MATRIX_CHECK_H_
#define PTOPE_UNIQUE_MATRIX_CHECK_H_

#include <armadillo>
#include <unordered_set>

#include "matrix_equiv.h"
#include "polytope_candidate.h"

namespace ptope {
template<class M>
class UMC {
	typedef std::unordered_set<M, ptope::MEquivHash, ptope::MEquivEqual>
		UniqueMSet;
public:
	bool
	operator()(const M & m);
	bool
	operator()(const PolytopeCandidate & p);
private:
	UniqueMSet _set;
};
typedef UMC<std::shared_ptr<const arma::mat>> UniqueMPtrCheck;
typedef UMC<arma::mat> UniqueMatrixCheck;
}
#endif

