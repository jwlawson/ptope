/*
 * matrix_hash.cc
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
#include "matrix_hash.h"

namespace ptope {
std::size_t
ColEquivSumHash::operator()(const arma::mat & m) const {
	std::size_t result{0};
	arma::uword n_col = m.n_cols;
	arma::uword n_row = m.n_rows;
	for(uint i = 0; i < n_col; ++i) {
		result += _vhash(n_row, m.colptr(i));
	}
	return result;
}
std::size_t
ColEquivSqSumHash::operator()(const arma::mat & m) const {
	std::size_t result{0};
	arma::uword n_col = m.n_cols;
	arma::uword n_row = m.n_rows;
	for(uint i = 0; i < n_col; ++i) {
		std::size_t h = _vhash(n_row, m.colptr(i));
		result += h * h;
	}
	return result;
}
std::size_t
ColEquivProdHash::operator()(const arma::mat & m) const {
	std::size_t result{0};
	arma::uword n_col = m.n_cols;
	arma::uword n_row = m.n_rows;
	for(uint i = 0; i < n_col; ++i) {
		result *= _vhash(n_row, m.colptr(i));
	}
	return result;
}
std::size_t
ColEquivDoublePlusHash::operator()(const arma::mat & m) const {
	static constexpr std::size_t R{2654435769};
	std::size_t result{0};
	arma::uword n_col = m.n_cols;
	arma::uword n_row = m.n_rows;
	for(uint i = 0; i < n_col; ++i) {
		result += (R + 2*_vhash(n_row, m.colptr(i)));
	}
	return result / 2;
}
std::size_t
VecHash::operator()(const arma::uword n_elem, const double * a) const {
	static constexpr std::size_t exp{31};
	std::size_t result{0};
	for(uint i = 0; i < n_elem; ++i) {
		result = result * exp + to_l(a[i]);
	}
	return result;
}
std::size_t
VecHash::operator()(const arma::vec & m) const {
	return operator()(m.size(), m.memptr());
}
}

