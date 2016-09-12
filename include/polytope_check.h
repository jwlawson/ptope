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

	struct Edge {
		vertex_index_t vertex;
		vector_index_t removed;
		Edge( vertex_index_t v, vector_index_t r ) : vertex(v), removed(r) {}
	};

	static constexpr vector_elem_t no_vertex = std::numeric_limits<vector_elem_t>::max();

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
	 * initial vertex, so this finds the other vertex along the edge.
	 *
	 * Such a vertex is not guaranteed to exist. In fact if no such vertex can be
	 * found then the provided gram matrix is not a hyperbolic polytope. In this
	 * case the function returns false.
	 *
	 * If a vertex is found, then it is copied into the provided vector and true
	 * is returned.
	 */
	vector_elem_t priv_edge_end( Edge const& edge, arma::mat const& gram,
			vector_t& vertex_out ) const;
	/**
	 * Find an initial elliptic subdiagram to use as initial vertex.
	 */
	void initial_vertex(PolytopeCandidate const& p, vector_t& output) const;
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
	 *
	 * The vector is constructed in the provided output vector. Ensure that the
	 * vector passed in has at least m_dimension elements. No size checking is
	 * done.
	 */
	void priv_vertex_from_edge(Edge const& cur_edge,
			vector_elem_t const vertex_index, vector_t& output) const;
	/**
	 * Check whether the given matrix is elliptic (i.e. positive definite).
	 */
	bool is_elliptic(arma::mat const& mat) const;
	/**
	 * Check whether the given matrix has a cholesky decomposition.
	 */
	bool priv_has_chol(arma::mat const& mat, arma::blas_int nrows,
		arma::blas_int ldmat) const;
	/**
	 * Copy the first num elements indexed by the indices vector from source_ptr
	 * to col_ptr.
	 */
	void priv_copy_submat_col( double * col_ptr, double const * source_ptr,
			vector_t const& indices, vector_index_t num ) const;
	/**
	 * Copy a submatrix from source to dest, only copying the elements indexed by
	 * the indices vector. Only the upper tirangular part of dest is written, but
	 * the matrix it points to is assumed to be a full (lddest x num_cols) matrix.
	 */
	void priv_copy_upper_triangle_submat( double * dest, double const * source,
			vector_t const& indices, vector_index_t num_cols, vector_index_t lddest,
			vector_index_t ldsource) const;
};

inline
void
PolytopeCheck::add_edges_from_vertex( vector_t const& vertex,
		vertex_index_t const vertex_ind, vector_elem_t const exclude) {
	for(vector_index_t i = 0, max = m_dimension; i < max; ++i) {
		if(vertex[i] == exclude) { continue; }
		_edge_queue.emplace( vertex_ind, i );
	}
}
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
	return priv_has_chol(mat, mat.n_rows, mat.n_cols);
}
// It might help to unroll this loop, but I'm not sure.
inline
void
PolytopeCheck::priv_copy_submat_col( double * col_ptr, double const * source_ptr,
			vector_t const& indices, vector_index_t num ) const {
	for( vector_index_t col_ind = 0; col_ind < num; ++col_ind ) {
		col_ptr[col_ind] = source_ptr[ indices(col_ind) ];
	}
}
inline
void
PolytopeCheck::priv_copy_upper_triangle_submat( double * dest,
		double const * source, vector_t const& indices, vector_index_t num_cols,
		vector_index_t lddest, vector_index_t ldsource ) const {
	for( vector_index_t col = 0; col < num_cols; ++col, dest += lddest ) {
		double const * next_source_col = source + indices( col ) * ldsource;
		priv_copy_submat_col( dest, next_source_col, indices, col + 1 );
	}
}
}
#endif

