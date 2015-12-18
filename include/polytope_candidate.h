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
/**
 * Candidate to be a polytope. Contains both the gram matrix and the vectors
 * which make up the polytope.
 *
 * Warning: This class is *not* thread safe for use in shared memory
 * applications.
 */
#pragma once
#ifndef PTOPE_POLYTOPE_CANDIDATE_H_
#define PTOPE_POLYTOPE_CANDIDATE_H_

#include <memory>

#include "vector_family.h"

namespace ptope {
class PolytopeCandidate {
	public:
		typedef std::shared_ptr<arma::mat> GramMatrix;
		/**
		 * Default constructor. No methods will work with an instance created using
		 * this. Just here for compatability.
		 */
		PolytopeCandidate();
		PolytopeCandidate(const PolytopeCandidate & p);
		/**
		 * Create a polytope candidate from an initial Gram matrix. The matrix can
		 * then be used to construct a set of vectors which form the normal vectors
		 * of the hyperplanes of the polytope. From there a change of corrdinates to
		 * the standard (orthonormal) basis is calculated.
		 */
		PolytopeCandidate(const GramMatrix & matrix);
		PolytopeCandidate(GramMatrix && matrix);
		PolytopeCandidate(const arma::mat & matrix);
		PolytopeCandidate(arma::mat && matrix);
		/**
		 * Create a polytope candidate from specified c-style arrays of specified
		 * size. Don't use unless you know what you are doing. The arrays are
		 * copied, and the PolytopeCandidate will not manage the pointers.
		 */
		PolytopeCandidate(const double * gram_ptr, int gram_size,
				const double * vector_ptr, int vector_dim, int no_vectors);
		PolytopeCandidate(std::initializer_list<std::initializer_list<double>> l);
		/**
		 * Given a vector of inner products with the basis vectors extend the
		 * polytope to include the new hyperplane defined by this vector.
		 */
		PolytopeCandidate
		extend_by_inner_products(const arma::vec & new_vector) const;
		/**
		 * Given a vector of inner products with the basis vectors extend the
		 * polytope to include the new hyperplane defined by this vector.
		 *
		 * The result is placed in the provided polytope candidate.
		 */
		bool
		extend_by_inner_products(PolytopeCandidate & result,
				const arma::vec & new_vector) const;
		/**
		 * Extend the polytope by a given normal vector.
		 */
		PolytopeCandidate
		extend_by_vector(const arma::vec & new_vector) const;
		/**
		 * Extend the polytope by a given normal vector, with the result in the
		 * provided polytope candidate.
		 */
		void
		extend_by_vector(PolytopeCandidate & result, const arma::vec & new_vector)
			const;
		/**
		 * Change the order of the vectors, so that a new basis consisting of those
		 * vectors specified by the given vector is used. This then constructs a new
		 * change of coordinates.
		 */
		void
		rebase_vectors(arma::uvec vec_indices);
		/**
		 * Return a copy of the polytope with vectors at indices a and b swapped.
		 */
		PolytopeCandidate
		swap_rebase(const arma::uword & a, const arma::uword & b) const;
		/**
		 * Get the signature of the polytope's gram matrix.
		 */
		std::pair<uint, uint>
		signature() const;
		/** Check whether the polytope is valid. */
		bool
		valid() const {
			return _valid;
		}
		std::size_t
		real_dimension() const;
		/**
		 * Return a reference to the polytope's gram matrix.
		 */
		const arma::mat &
		gram() const {
			return *_gram;
		}
		std::shared_ptr<const arma::mat>
		gram_ptr() const {
			return static_cast<std::shared_ptr<const arma::mat>>(_gram);
		}
		/**
		 * Get a reference to the polytope's vector family.
		 */
		const VectorFamily &
		vector_family() const {
			return _vectors;
		}
		/**
		 * Save the polytope candidate to stream.
		 */
		void
		save(std::ostream & os);
		/**
		 * Read the polytope from stream.
		 */
		void
		load(std::istream & is);
		/**
		 * Swap the resources of this polytope with the one provided.
		 */
		void
		swap(PolytopeCandidate & p);
		PolytopeCandidate &
		operator=(const PolytopeCandidate & p);
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

		/**
		 * Calculate the vector which gives the specified inner products.
		 * Returns true if the vector is valid, false if it is invalid.
		 */
		bool
		vector_from_inner_products(const arma::vec & inner_vector) const;
};
}
#endif

