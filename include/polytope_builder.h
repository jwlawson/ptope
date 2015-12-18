/*
 * polytope_builder.h
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
/*
 * Tries to build a compact polytope out of a gram matrix which is not a
 * (compact) polytope. PolytopeCheck looks for edges which do not have two
 * elliptic vertices, so the idea is to add additional hyperplanes to the
 * polytope to close off any unending edges.
 */
#pragma once
#ifndef PTOPE_POLYTOPE_BUILDER_H_
#define PTOPE_POLYTOPE_BUILDER_H_

#include "inner_product_vectors.h"
#include "polytope_candidate.h"
#include "polytope_check.h"

namespace ptope {
class PolytopeBuilder {
	static constexpr int max_depth = 4;
	typedef std::vector<PolytopeCandidate> Results;
	typedef std::vector<InnerProductVectors> VectorSet;
public:
	typedef Results::iterator iterator;
	typedef Results::const_iterator const_iterator;

	PolytopeBuilder(int size)
		: _vectors(max_depth, InnerProductVectors(size)) {}
	/**
	 * Add vectors to the supplied polytope to construct polytopes.
	 */
	iterator
	build(const PolytopeCandidate & p);
	/**
	 * Get an iterator to the end of the results.
	 */
	iterator
	end();
private:
	void
	compute_polytopes(const PolytopeCandidate & p);
	/**
	 * Find inner product vectors which make an unending edge elliptic.
	 */
	void
	find_elliptic_vector();
	/**
	 * Check if an inner_product vector is compatible with the gram matrix.
	 */
	void
	compatible_vector();
	void
	problem_edge();
	void
	extend_elliptic_vector();
	void
	extend_till_polytope(const PolytopeCandidate & p, int depth);

	/** Stores the resulting polytopes. */
	PolytopeCheck _chk;
	Results _results;
	VectorSet _vectors;
	arma::vec _ext_vec;
};
}
#endif

