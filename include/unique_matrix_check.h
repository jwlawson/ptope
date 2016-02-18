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

#include "boost/bloom_filter/basic_bloom_filter.hpp"

#include "matrix_equiv.h"
#include "matrix_hash.h"
#include "polytope_candidate.h"

namespace ptope {
class UniquePCCheck {
	typedef std::unordered_set<arma::mat, ptope::MEquivHash, ptope::MEquivEqual>
		UniqueMSet;
public:
	bool
	operator()(const arma::mat & m);
	bool
	operator()(const PolytopeCandidate & p);
private:
	UniqueMSet _set;
};
class BloomPCCheck {
static constexpr std::size_t max_size = std::numeric_limits<std::size_t>::max();
static constexpr std::size_t default_size = 8589934592ull;
typedef boost::mpl::vector< ColEquivSumHash,
														ColEquivProdHash,
														ColEquivSqSumHash,
														ColEquivDoublePlusHash > HashVector;
typedef boost::bloom_filters::basic_bloom_filter<arma::mat, default_size,
				HashVector> Filter;
public:
	bool
	operator()(const arma::mat & m);
	bool
	operator()(const PolytopeCandidate & p);
private:
	Filter _filter;
};
}
#endif

