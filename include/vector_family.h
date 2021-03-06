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
		 * Construct a vector family from the given c-style array, for the given
		 * number of vector_dim dimensional vectors.
		 */
		VectorFamily(const double * vector_ptr, int vector_dim, int no_vectors);
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
		 * Copy the provided vector family, and add an additional vector.
		 */
		void
		copy_and_add_vector(const VectorFamily & vf, const arma::vec & vec);
		/**
		 * Copy the provided vector family, and add an additional vector.
		 */
		void
		copy_and_add_first_hyperbolic_vector(const VectorFamily & vf,
				const arma::vec & vec);
		/**
		 * Swap two vectors in the family.
		 */
		void
		swap(const arma::uword a, const arma::uword b) {
			_vectors.swap_cols(a, b);
		}
		/**
		 * Get a vector from the family. The unsafe version gets a vector which
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
		 * Get a pointer to a vector in the family.
		 */
		double const *
		get_ptr(const arma::uword index) const;
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
		/**
		 * Get the first columns which make up a basis of the vector family.
		 *
		 * In fact you get one fewer vector than a basis, as the requirement that
		 * the vectors be unit gives the final coordinate.
		 */
		const arma::mat
		first_basis_cols() const;
		/**
		 * Save the vector family to the stream.
		 */
		void
		save(std::ostream & os, arma::file_type type) const;
		/**
		 * Load a saved vector family from stream.
		 */
		void
		load(std::istream & os, arma::file_type type);
		/**
		 * Swap resources with provided VectorFamily.
		 */
		void
		swap(VectorFamily & vf);
	private:
		arma::mat _vectors;
};
inline
arma::vec
VectorFamily::unsafe_get(const arma::uword index) const {
	return _vectors.unsafe_col(index);
}
inline
arma::vec
VectorFamily::get(const arma::uword index) const {
	return _vectors.col(index);
}
inline
double const *
VectorFamily::get_ptr(const arma::uword index) const {
	return _vectors.colptr(index);
}
inline
const arma::mat
VectorFamily::first_basis_cols() const  {
	arma::uword size = dimension() - 1;
	return _vectors.head_cols(size);
}
inline
void
VectorFamily::save(std::ostream & os, arma::file_type type) const {
	_vectors.save(os, type);
}
inline
void
VectorFamily::load(std::istream & is, arma::file_type type) {
	_vectors.load(is, type);
}
inline
void
VectorFamily::swap(VectorFamily & vf) {
	_vectors.swap(vf._vectors);
}
}
#endif

