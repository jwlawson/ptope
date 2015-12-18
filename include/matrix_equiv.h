/*
 * matrix_equiv.h
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
#ifndef PTOPE_MATRIX_H_
#define PTOPE_MATRIX_H_

#include <armadillo>
#include <memory>

namespace ptope {
struct MEquivHash {
	std::size_t
	operator()(const arma::mat & m) const;
	std::size_t
	operator()(const std::shared_ptr<arma::mat> & p) const {
		return operator()(*p);
	}
};
struct MEquivEqual {
	/**
	 * Check whether the two provided matrices are the same up to some
	 * permutation.
	 */
	bool
	operator()(const arma::mat & lhs, const arma::mat & rhs) const;
	bool
	operator()(const std::shared_ptr<arma::mat> & lhs,
			const std::shared_ptr<arma::mat> & rhs) const {
		return operator()(*lhs, *rhs);
	}
private:
	/**
	 * Get the column sums for the provided matrix.
	 */
	void
	sums(const arma::mat & m, std::vector<double> & result) const;
	/**
	 * Find which columns of the matrix can possibly be mapped to which columns in
	 * the other matrix. Thislimits the possible permutations to check.
	 */
	void
	compatible_vectors(const arma::mat & lhs, const arma::mat & rhs,
			std::vector<std::vector<arma::uword>> & result) const;
	/**
	 * Recursively construct a possible permutation of the vectors and once done
	 * check whether that permutation takes one vector to the other.
	 */
	bool
	check_permutation(const arma::mat & lhs, const arma::mat & rhs,
			std::vector<arma::uword> & perm, std::size_t index,
			const std::vector<std::vector<arma::uword>> & comp) const;
};
}
#endif

