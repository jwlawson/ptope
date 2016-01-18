/*
 * polytope_rebaser.h
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
#ifndef PTOPE_POLYTOPE_REBASER_H_
#define PTOPE_POLYTOPE_REBASER_H_

#include "polytope_candidate.h"

namespace ptope {
class PolytopeRebaser {
class PermIter {
public:
	PermIter(int size, int max);
	bool
	has_next();
	const arma::uvec &
	next();
private:
	arma::uvec _next;
	std::vector<arma::uword> _progress;
	const arma::uword _progress_max;
	const arma::uword _diff;
	bool _has_next;

	void
	increment_progress();
	void
	increment(const std::size_t & ind);
	void
	compute_next();
};
public:
	PolytopeRebaser(const PolytopeCandidate & p);
	/**
	 * Return true if a subsequent call to next() will return a valid polytope
	 * candidate.
	 */
	bool
	has_next();
	/**
	 * Return the next rebased polytope.
	 */
	const PolytopeCandidate &
	next();
private:
	PolytopeCandidate _initial;
	PermIter _perm;

};
}
#endif

