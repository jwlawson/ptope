/*
 * parabolic_check.h
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
 * Functor to check whether a given connected symmetric matrix is parabolic.
 *
 * A connected matrix is parabolic if it is positive semi-definite but not
 * positive definite. Equivalently all eigenvalues are non-negative, with at
 * least one being 0.
 */
#pragma once
#ifndef PTOPE_PARABOLIC_CHECK_H_
#define PTOPE_PARABOLIC_CHECK_H_

#include "polytope_candidate.h"

namespace ptope {
class ParabolicCheck {
	public:
		ParabolicCheck() {}
		/**
		 * Check whether the given symmetric matrix is parabolic.
		 * Returns true if it is parabolic.
		 */
		bool operator()(const arma::mat & m);
		bool operator()(const PolytopeCandidate & p);
	private:
		arma::vec evalues;
};
}
#endif

