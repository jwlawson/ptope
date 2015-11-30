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
	PolytopeCandidate
	next();
private:
	PolytopeCandidate _initial;
	arma::uword _basis_ind;
	const arma::uword _basis_max;
	arma::uword _swap_ind;
	const arma::uword _swap_max;
};
}
#endif

