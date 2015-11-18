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
void VectorFamily::add_vector(const arma::mat & vec) {
	arma::uword last_col = _vectors.n_cols;
	_vectors.set_size(_vectors.n_rows, _vectors.n_cols + 1);
	for(arma::uword i = 0; i < _vectors.n_rows; ++i) {
		_vectors(i, last_col) = vec(i);
	}
}
arma::vec VectorFamily::unsafe_get(const arma::uword index) const {
	return _vectors.unsafe_col(index);
}
arma::vec VectorFamily::get(const arma::uword index) const {
	return _vectors.col(index);
}
}

