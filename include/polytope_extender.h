/*
 * polytope_extender.h
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
/**
 * Class to extend a given polytope in every possible way, given certain angles
 * to extend the polytope with.
 *
 * Uses a java iterator style interface: has_next() to check if there are
 * further polytopes and next() to generate and return the next one.
 */
#pragma once
#ifndef PTOPE_POLYTOPE_EXTENDER_H_
#define PTOPE_POLYTOPE_EXTENDER_H_

#include "polytope_candidate.h"

namespace ptope {
class PolytopeExtender {
	static std::vector<double> __default_inner;
	public:
		PolytopeExtender(const PolytopeCandidate & initial_polytope);
		PolytopeExtender(PolytopeCandidate && initial_polytope);
		PolytopeExtender(const PolytopeCandidate & initial_polytope,
				const std::vector<uint> & angles);
		PolytopeExtender(PolytopeCandidate && initial_polytope,
				const std::vector<uint> & angles);
		/**
		 * Check whether a subsequent call to next() will return a valid polytope.
		 */
		bool
		has_next();
		/**
		 * Compute and return the next extended polytope.
		 */
		const PolytopeCandidate &
		next();
		/**
		 * Set the default angles used to extend the polytopes.
		 */
		static void
		set_default_angles(const std::vector<uint> & angles);
	private:
		PolytopeCandidate _initial;
		std::vector<double> _inner_products;
		arma::vec _extend_by;
		std::vector<std::size_t> _progress;
		std::size_t _progress_max;
		bool _has_next;
		PolytopeCandidate _next;
		bool _computed_next;

		/**
		 * Increment the _progress counter, so next call to next() will provide the
		 * next polytope.
		 */
		void
		increment_progress();
		bool
		compute_next();
};
}
#endif

