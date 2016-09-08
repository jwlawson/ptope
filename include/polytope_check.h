/*
 * polytope_check.h
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
#ifndef PTOPE_POLYTOPE_CHECK_H_
#define PTOPE_POLYTOPE_CHECK_H_

#include "comparator.h"
#include "polytope_candidate.h"
#include "vector_set.h"

#include <algorithm>
#include <queue>
#include <set>
#include <vector>

#include "boost/pool/pool_alloc.hpp"

namespace ptope {
class PolytopeCheck {
private:
public:
	PolytopeCheck();
	/**
	 * Check whether a given polytope candidate is actually a compact polytope.
	 * Returns true if it is a compact polytope.
	 */
	bool operator()(PolytopeCandidate const& p);
private:
	typedef uint_fast16_t vector_elem_t;
	typedef arma::Col<vector_elem_t> vector_t;
	typedef std::size_t vertex_index_t;
	typedef uint_fast16_t vector_index_t;

	/** Value to return when no vertex can be found. */
	static constexpr vector_elem_t
	no_vertex = std::numeric_limits<vector_elem_t>::max();

	struct Edge {
		vertex_index_t vertex;
		vector_index_t removed;
	};


	typedef boost::pool_allocator<Edge,
					boost::default_user_allocator_new_delete,
					boost::details::pool::null_mutex> EdgeAllocator;
	typedef std::deque<Edge, EdgeAllocator> EdgeContainer;
	typedef std::queue<Edge, EdgeContainer> EdgeQueue;

	VectorSet<vector_elem_t> _visited_vertices;
	EdgeQueue _edge_queue;
	vector_index_t m_dimension;
	/**
	 * Find the vertex at the end of an edge. Each edge is constructed from an
	 * initial vertex, so this finds the other vertex along the edge. If no vertex
	 * is found then PolytopeCheck::no_vertex is returned.
	 *
	 * The index returned is the additional vertex index needed to complete the
	 * edge to a vertex.
	 */
	vector_elem_t
	find_edge_end(Edge const& edge, arma::mat const& p) const;
	/**
	 * Find an initial elliptic subdiagram to use as initial vertex.
	 */
	vector_t initial_vertex(PolytopeCandidate const& p) const;
	/**
	 * Use the provided vertes to add all edges adjacent to that vertex to the
	 * queue, except for the edge specified by the excluded index. The excluded
	 * index corresponds to the edge which first visited this vertex.
	 */
	void
	add_edges_from_vertex(vector_t const& vertex, vertex_index_t const vertex_ind,
			vector_elem_t const exclude);
	/**
	 * Construct the vector representing the vertex which lies at the end of the
	 * specified edge with the given additional index.
	 */
	vector_t get_vertex_from_edge(Edge const& cur_edge,
			vector_elem_t const vertex_index) const;
	/**
	 * Check whether the given matrix is elliptic (i.e. positive definite).
	 */
	bool is_elliptic(arma::mat const& mat) const;
	/**
	 * Check whether the given matrix has a cholesky decomposition.
	 */
	bool has_chol(arma::mat const& mat) const;
	/**
	 * Get the vertex at the specified index.
	 *
	 * The VectorSet does not have a const at(index) method, as it cannot
	 * guarantee that the provided vector will not be altered.
	 *
	 * !Any use of the returned vector should be const!
	 */
	vector_t priv_const_vertex_at( vertex_index_t const index ) const;
};
/*
 * A matrix is positive definite iff all its eigen values are positive. However
 * finding the eigenvalues of a matrix is computationally hard. Equivalently
 * there are two other definitions which help:
 *
 *  -Sylvesters criterion: A matrix is positive definite iff all determinants of
 *    leading minors are positive. This could be checked with Gaussian
 *    elimination to get the matrix in upper triangular form (using BLAS/LAPACK
 *    LU decomposition) so the determinants can be read off the diagonal entries.
 *
 *  -Cholesky decomp: A matrix is positive definite iff it has a Cholesky
 *  	decomposition. This can easily be checked by trying to compute a Cholesky
 *  	decomp (using BLAS/LAPACK) and seeing if it successful.
 */
inline
bool
PolytopeCheck::is_elliptic(arma::mat const& mat) const {
	bool success = has_chol(mat);
	return success;
}
inline
typename PolytopeCheck::vector_t
PolytopeCheck::priv_const_vertex_at( vertex_index_t const index ) const {
	vector_elem_t const * cvptr = _visited_vertices.ptr_at( index );
	vector_elem_t * ptr = const_cast<vector_elem_t *>( cvptr );
	return vector_t( ptr, m_dimension, false, true );
}
}
#endif

