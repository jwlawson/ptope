/*
 * vector_family.cc
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
#include "vector_family.h"

namespace ptope {
VectorFamily::VectorFamily(const arma::mat & mat)
	: _vectors(mat) {}
VectorFamily::VectorFamily(arma::mat && mat)
	: _vectors(mat) {}
VectorFamily::VectorFamily(const double * vector_ptr, int vector_dim,
		int no_vectors)
	: _vectors(vector_ptr, vector_dim, no_vectors) {}

void
VectorFamily::add_vector(const arma::mat & vec) {
	arma::uword last_col = _vectors.n_cols;
	_vectors.resize(_vectors.n_rows, _vectors.n_cols + 1);
	for(arma::uword i = 0, max = _vectors.n_rows; i < max; ++i) {
		_vectors(i, last_col) = vec(i);
	}
}
void
VectorFamily::add_first_hyperbolic_vector(const arma::mat & vec) {
	arma::uword last_col = _vectors.n_cols;
	arma::uword last_row = _vectors.n_rows;
	_vectors.resize(_vectors.n_rows + 1, _vectors.n_cols + 1);
	arma::uword size = vec.size();
	for(arma::uword i = 0; i < size; ++i) {
		_vectors(last_row, i) = 0;
		_vectors(i, last_col) = vec(i);
	}
}
void
VectorFamily::copy_and_add_vector(const VectorFamily & vf,
		const arma::vec & vec) {
	const arma::uword & last_col = vf._vectors.n_cols;
	const arma::uword & last_row = vf._vectors.n_rows;
	_vectors.set_size(last_row, last_col + 1);
	_vectors.head_cols(last_col) = vf._vectors;
	for(arma::uword i = 0, max = last_row; i < max; ++i) {
		_vectors(i, last_col) = vec(i);
	}
}
void
VectorFamily::copy_and_add_first_hyperbolic_vector(const VectorFamily & vf,
		const arma::vec & vec) {
	const arma::uword & last_col = vf._vectors.n_cols;
	const arma::uword & last_row = vf._vectors.n_rows;
	_vectors.set_size(last_row + 1, last_col + 1);
	_vectors.submat(0, 0, last_row - 1, last_col - 1) = vf._vectors;
	const arma::uword & size = vec.size();
	for(arma::uword i = 0; i < size; ++i) {
		_vectors(last_row, i) = 0;
		_vectors(i, last_col) = vec(i);
	}
}
arma::vec
VectorFamily::unsafe_get(const arma::uword index) const {
	return _vectors.unsafe_col(index);
}
arma::vec
VectorFamily::get(const arma::uword index) const {
	return _vectors.col(index);
}
const arma::mat
VectorFamily::first_basis_cols() const  {
	arma::uword size = dimension() - 1;
	return _vectors.head_cols(size);
}
void
VectorFamily::save(std::ostream & os, arma::file_type type) {
	_vectors.save(os, type);
}
void
VectorFamily::load(std::istream & is, arma::file_type type) {
	_vectors.load(is, type);
}
}


