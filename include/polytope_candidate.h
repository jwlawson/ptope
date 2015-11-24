/*
 * polytope_candidate.h
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
#ifndef PTOPE_POLYTOPE_CANDIDATE_H_
#define PTOPE_POLYTOPE_CANDIDATE_H_

#include "vector_family.h"

namespace ptope {
class PolytopeCandidate {
	public:
		typedef arma::mat GramMatrix;
		/**
		 * Create a polytope candidate from an initial Gram matrix. The matrix can
		 * then be used to construct a set of vectors which form the normal vectors
		 * of the hyperplanes of the polytope. From there a change of corrdinates to
		 * the standard (orthonormal) basis is calculated.
		 */
		PolytopeCandidate(const GramMatrix & matrix);
		PolytopeCandidate(GramMatrix && matrix);
		/**
		 * Given a vector of inner products with the basis vectors extend the
		 * polytope to include the new hyperplane defined by this vector.
		 */
		PolytopeCandidate
		extend_by_inner_products(const arma::vec & new_vector);
		/**
		 * Change the order of the vectors, so that a new basis consisting of those
		 * vectors specified by the given vector is used. This then constructs a new
		 * change of coordinates.
		 */
		void
		rebase_vectors(arma::uvec vec_indices);
		/**
		 * Get the signature of the polytope's gram matrix.
		 */
		std::pair<uint, uint>
		signature() const;
		/** Check whether the polytope is valid. */
		bool
		valid() {
			return _valid;
		}
		/**
		 * Print matrix and vectors to output stream.
		 */
		friend
		std::ostream &
		operator<<(std::ostream & os, const PolytopeCandidate & poly);
	private:
		static PolytopeCandidate InValid;
		/** Gram matrix of the polytope. */
		GramMatrix _gram;
		/** Set of normal vectors of hyperplanes.
		 *
		 * The first _vectors.dimension() vectors are used as a basis of the vector
		 * space.
		 */
		VectorFamily _vectors;
		/**
			* Transpose of matrix of basis vectors.
			*
			* If R is the matrix of (d+1)-dim vectors, then this is R^T.
			* This is useful as the construction of a vector A given the inner_product
			* vector B is the same as solving (R^T).A = B.
			*/
		arma::mat _basis_vecs_trans;
		/** Is the set of vectors hyperbolic or still real? */
		bool _hyperbolic;
		/** Check whether the polytope is valid */
		bool _valid;
		/** Private default constructor to allow creation of Invalid. Don't use! */
		PolytopeCandidate() : _vectors(arma::mat()), _valid(false) {}
};
}
#endif

