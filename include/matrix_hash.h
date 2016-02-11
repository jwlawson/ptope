/*
 * matrix_hash.h
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
#ifndef PTOPE_MATRIX_HASH_H_
#define PTOPE_MATRIX_HASH_H_

#include <armadillo>

namespace ptope {
struct VecHash {
	std::size_t
	operator()(const arma::uword n_elem, const double * a) const;
	std::size_t
	operator()(const arma::vec & m) const;
private:
	long
	inline
	to_l(const double & d) const {
		return std::lround(d * 1e5);
	}
};
struct ColEquivSumHash {
	std::size_t
	operator()(const arma::mat & m) const;
private:
	VecHash _vhash;
};
struct ColEquivSqSumHash {
	std::size_t
	operator()(const arma::mat & m) const;
private:
	VecHash _vhash;
};
struct ColEquivProdHash {
	std::size_t
	operator()(const arma::mat & m) const;
private:
	VecHash _vhash;
};
struct ColEquivDoublePlusHash {
	std::size_t
	operator()(const arma::mat & m) const;
private:
	VecHash _vhash;
};
}
#endif

