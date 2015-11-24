/*
 * vector_family.h
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
#ifndef PTOPE_VECTOR_FAMILY_H_
#define PTOPE_VECTOR_FAMILY_H_

#include <armadillo>

namespace ptope {
class VectorFamily {
	public:
		/**
		 * Construct a vector family from a matrix. Each column of the matrix
		 * corresponds to a vector in the family. The initial matrix specifies d
		 * real vectors of dimension d.
		 */
		VectorFamily(const arma::mat & matrix);
		VectorFamily(arma::mat && matrix);
		/**
		 * Add a vector the family. The vector is assumed to be the same size as all
		 * others in the family.
		 */
		void
		add_vector(const arma::mat & vec);
		/**
		 * Add a vector the family. The vector is assumed to be one entry smaller
		 * than all others in the family, with the final entry given separately as
		 * im_coord.
		 */
		void
		add_first_hyperbolic_vector(const arma::mat & vec);
		/**
		 * Swap two vectors in the family.
		 */
		void
		swap(const arma::uword a, const arma::uword b) {
			_vectors.swap_cols(a, b);
		}
		/**
		 * Get a vetor from the family. The unsafe version gets a vector which
		 * points to the memory of the vectors in the vector family. Only use if you
		 * know what you are doing.
		 */
		arma::vec
		unsafe_get(const arma::uword index) const;
		/**
		 * Get a vector from the family. The safe version returns a copy of the
		 * vector and so is safe to change.
		 */
		arma::vec
		get(const arma::uword index) const;
		/**
		 * Get the dimension of each vector in the family.
		 */
		arma::uword
		dimension() const {
			return _vectors.n_rows;
		}
		/**
		 * Get the number of vectors in the family.
		 */
		arma::uword
		size() const {
			return _vectors.n_cols;
		}
		/**
		 * Get a const reference to the underlying matrix used to store the
		 * vectors in the family. 
		 */
		const arma::mat &
		underlying_matrix() const {
			return _vectors;
		}
		const arma::mat
		first_basis_cols() const;
	private:
		arma::mat _vectors;
};
}
#endif

