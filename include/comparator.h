/*
 * calc.h
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
 * Provides calculations and simple functions which varioous classes need to
 * evaluate. For example min_cos_angle and mik_inner_product.
 */
#pragma once
#ifndef _PTOPE_COMPARATOR_H_
#define _PTOPE_COMPARATOR_H_

#include <armadillo>

namespace ptope {
namespace comparator {
/** Error tolerance in double calculations. */
constexpr double error = 1e-10;
/** Comparator for doubles with tolerance of calc::error. */
struct DoubleLess {
bool
operator()(double const & lhs, double const & rhs) const {
	return (lhs + error < rhs);
}
};
/** Comparator for arma::uvec using lexicographic comparison. */
struct UVecLess {
bool
operator()(arma::uvec const & lhs, arma::uvec const & rhs) const {
	return std::lexicographical_compare(lhs.begin(), lhs.end(),
			rhs.begin(), rhs.end());
}
bool
operator()(arma::uvec const * const lhs, arma::uvec const * const rhs) const {
	return operator()(*lhs, *rhs);
}
};
/** Comparator for arma::vec using lexicographic comparison. */
struct VecLess {
bool 
operator()(arma::vec const & lhs, arma::vec const & rhs) const {
	return std::lexicographical_compare(lhs.begin(), lhs.end(),
			rhs.begin(), rhs.end(), _dless);
}
bool
operator()(arma::vec const * const lhs, arma::vec const * const rhs) const {
	return operator()(*lhs, *rhs);
}
private:
DoubleLess _dless;
};

}
}
#endif

