/*
 * matrix_equiv.cc
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
#include "matrix_equiv.h"

namespace ptope {
std::size_t
MEquivHash::operator()(const arma::mat & m) const {
	std::size_t result = 0;
	std::size_t R = 1779033703;
	for(arma::uword i = 0, max = m.n_cols; i < max; ++i) {
		std::size_t c_sum = 0;
		std::size_t c_sq_sum = 0;
		for(auto cit = m.begin_col(i), end = m.end_col(i); cit != end; ++cit) {
			std::size_t val = std::lround(*cit * 10000);
			c_sum += R + 2*val;
			c_sq_sum += val * val;
		}
		result += c_sum * c_sum + (R + 2*c_sq_sum);
	}
	result /= 2;
	return result;
}
bool
MEquivEqual::operator()(const arma::mat & lhs, const arma::mat & rhs) const {
	bool result = false;
	if(lhs.n_cols == rhs.n_cols) {
		__perm.reserve(lhs.n_cols);
		compatible_vectors(lhs, rhs, __comp);
		result = check_permutation(lhs, rhs, __perm, 0, __comp);
		__perm.clear();
	}
	return result;
}
void
MEquivEqual::sums(const arma::mat & m, std::vector<double> & out) const {
	out.reserve(m.n_cols);
	for(arma::uword i = 0, max = m.n_cols; i < max; ++i) {
		double c_sum = 0;
		for(auto cit = m.begin_col(i), end = m.end_col(i); cit != end; ++cit) {
			c_sum += *cit;
		}
		out[i] = c_sum;
	}
}
void
MEquivEqual::compatible_vectors(const arma::mat & lhs, const arma::mat & rhs,
		std::vector<std::vector<arma::uword>> & result) const {
	for(auto & vec : result) {
		vec.clear();
	}
	for(std::size_t i = result.size(), max = lhs.n_cols; i < max; ++i) {
		result.push_back(std::vector<arma::uword>());
	}
	sums(lhs, __l_sum);
	sums(rhs, __r_sum);
	for(arma::uword j = 0, max = lhs.n_cols; j < max; ++j) {
		for(arma::uword i = 0; i < max; ++i) {
			if(_d_eq(__l_sum[j], __r_sum[i])) result[j].push_back(i);
		}
	}
}
bool
MEquivEqual::check_permutation(const arma::mat & lhs, const arma::mat & rhs,
		std::vector<arma::uword> & perm, std::size_t index,
		const std::vector<std::vector<arma::uword>> & comp) const {
	bool result = false;
	if(index == lhs.n_cols) {
		// Have complete permutation
		__p_cols.set_size(rhs.n_rows, rhs.n_cols);
		for(arma::uword i = 0, max = lhs.n_cols; i < max; ++i) {
			__p_cols.col(i) = rhs.col(perm[i]);
		}
		__permuted.set_size(rhs.n_rows, rhs.n_cols);
		for(arma::uword i = 0, max = lhs.n_cols; i < max; ++i) {
			__permuted.row(i) = __p_cols.row(perm[i]);
		}
		bool equal = true;
		for(arma::uword i = 0, max = __permuted.size(); i < max && equal; ++i) {
			equal = _d_eq(__permuted(i), lhs(i));
		}
		result = equal;
	} else {
		// Add another entry to the permutation
		for(const arma::uword & map : comp[index]) {
			if(result) break;
			if(std::find(perm.begin(), perm.begin() + index, map)
					== perm.begin() + index) {
				perm[index] = map;
				result = check_permutation(lhs, rhs, perm, index + 1, comp);
			}
		}
	}
	return result;
}
bool
MColPermEquiv::operator()(const arma::mat & lhs, const arma::mat & rhs) const {
	bool result = false;
	if(lhs.n_cols == rhs.n_cols) {
		__perm.reserve(lhs.n_cols);
		compatible_vectors(lhs, rhs, __comp);
		result = check_permutation(lhs, rhs, __perm, 0, __comp);
		__perm.clear();
	}
	return result;
}
void
MColPermEquiv::sums(const arma::mat & m, std::vector<double> & out) const {
	out.reserve(m.n_cols);
	for(arma::uword i = 0, max = m.n_cols; i < max; ++i) {
		double c_sum = 0;
		for(auto cit = m.begin_col(i), end = m.end_col(i); cit != end; ++cit) {
			c_sum += *cit;
		}
		out[i] = c_sum;
	}
}
void
MColPermEquiv::compatible_vectors(const arma::mat & lhs, const arma::mat & rhs,
		std::vector<std::vector<arma::uword>> & result) const {
	for(auto & vec : result) {
		vec.clear();
	}
	for(std::size_t i = result.size(), max = lhs.n_cols; i < max; ++i) {
		result.push_back(std::vector<arma::uword>());
	}
	sums(lhs, __l_sum);
	sums(rhs, __r_sum);
	for(arma::uword j = 0, max = lhs.n_cols; j < max; ++j) {
		for(arma::uword i = 0; i < max; ++i) {
			if(_d_eq(__l_sum[j], __r_sum[i])) result[j].push_back(i);
		}
	}
}
bool
MColPermEquiv::check_permutation(const arma::mat & lhs, const arma::mat & rhs,
		std::vector<arma::uword> & perm, std::size_t index,
		const std::vector<std::vector<arma::uword>> & comp) const {
	bool result = false;
	if(index == lhs.n_cols) {
		// Have complete permutation
		__permuted.set_size(rhs.n_rows, rhs.n_cols);
		for(arma::uword i = 0, max = lhs.n_cols; i < max; ++i) {
			__permuted.col(i) = rhs.col(perm[i]);
		}
		bool equal = true;
		for(arma::uword i = 0, max = __permuted.size(); i < max && equal; ++i) {
			equal = _d_eq(__permuted(i), lhs(i));
		}
		result = equal;
	} else {
		// Add another entry to the permutation
		for(const arma::uword & map : comp[index]) {
			if(result) break;
			if(std::find(perm.begin(), perm.begin() + index, map)
					== perm.begin() + index) {
				perm[index] = map;
				result = check_permutation(lhs, rhs, perm, index + 1, comp);
			}
		}
	}
	return result;
}
}

