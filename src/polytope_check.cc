/*
 * polytope_check.cc
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
#include "polytope_check.h"

namespace ptope {
namespace {
arma::mat __elliptic_check_tmp;
arma::uvec __indices_cached;
arma::mat __vertex_submat;
}
PolytopeCheck::PolytopeCheck()
	: _last_edge(0, arma::uvec()) {}
/*
 * The Gram matrix of a polytope contains all the information to determine
 * whether or not it is compact - which is really what this method is checking.
 * The basic idea is to start with one vertex and consider all sides at that
 * vertex. If the polytope is compact, then each of those sides will have a
 * second vertex somewhere along them. Do a breadth-first search along all such
 * edges from all vertices to determine whether there is any edge which does not
 * have a second vertex along it.
 *
 * It is a result of Vinberg that a k-dim face in the polytope corresponds to a
 * (d-k)x(d-k) submatrix of the gram matrix. Hence vertices are dxd submatrices
 * and edges are (d-1)x(d-1) submatrices. It is also true that for compact
 * polytopes each vertex must be elliptic, as the vertex is contained within
 * hyperbolic space, so the reflection group generated by reflections through
 * the planes intersecting at the vertex is finite. (If the vertex were ideal
 * the group is parabolic).
 *
 * Therefore this check starts with an elliptic dxd submatrix of the gram
 * matrix and considers all possible (d-1)x(d-1) submatrices of that submatrix.
 * These are the edges, so adding some other row/column of the gram matrix to
 * get a dxd matrix will give something that might possibly be a vertex. To
 * check whether it is or not it suffices to see if that matrix is elliptic. If
 * so then we have the other vertex along the side, otherwise try a different
 * row/column.
 *
 * Eventually either all edges have been checked and found to contain two
 * vertices so the matrix is a polytope, or one edge will be found to only
 * contain one vertex so the matrix is not a polytope.
 */
bool
PolytopeCheck::operator()(const PolytopeCandidate & p) {
	_last_polytope = &p;

	const arma::uvec init_v = initial_vertex(p);
	for(arma::uword i = 0, max = init_v.size(); i < max; ++i) {
		_edge_queue.emplace(init_v(i), construct_edge(init_v, i));
	}
	_visited_vertices.emplace(std::move(init_v));
	bool result = true;
	while(result && !_edge_queue.empty()) {
		_last_edge = std::move(_edge_queue.front());
		_edge_queue.pop();
		if(!handle_edge(_last_edge, p) ) result = false;
	}
	/* _last_edge uses memory obtained from the pool. We are about to purge the
	 * pool, so need to reconstruct _last_edge in different memory. */
	Edge tmp(_last_edge);
	_last_edge = std::move(tmp);
	_visited_vertices.clear();
	while(!_edge_queue.empty()) _edge_queue.pop();
	boost::singleton_pool<boost::fast_pool_allocator_tag, sizeof(Edge)>::purge_memory();
	boost::singleton_pool<boost::fast_pool_allocator_tag, sizeof(arma::uvec)>::purge_memory();
	return result;
}
bool
PolytopeCheck::resume(const PolytopeCandidate & p) {
	_last_polytope = &p;

	if(!handle_edge(_last_edge, p) ) return false;

	while(!_edge_queue.empty()) {
		_last_edge = std::move(_edge_queue.front());
		_edge_queue.pop();
		if(!handle_edge(_last_edge, p) ) return false;
	}
	return false;
}
const arma::subview_elem2<double, arma::Mat<unsigned long long>,
			arma::Mat<unsigned long long> >
PolytopeCheck::last_edge() const {
	return _last_polytope->gram().submat(_last_edge.edge, _last_edge.edge);
}
bool
PolytopeCheck::handle_edge(const Edge & e, const PolytopeCandidate & p) {
	const arma::uword next_vert_ind = find_edge_end(e, p);
	if(next_vert_ind == no_vertex) {
		return false;
	}
	const arma::uvec new_vert = get_vertex_from_edge(e, next_vert_ind);
	if(vertex_unvisited(new_vert)) {
		add_edges_from_vertex(new_vert, next_vert_ind);
		_visited_vertices.emplace(std::move(new_vert));
	}
	return true;
}
arma::uword
PolytopeCheck::find_edge_end(const Edge & edge,
		const PolytopeCandidate & p) const {
	const arma::mat & gram = p.gram();
	const arma::uword & last_entry = edge.edge.size();
	const arma::uword & size = last_entry + 1;
	__indices_cached.set_size(size);
	__vertex_submat.set_size(size, size);
	for(arma::uword k = 0; k < last_entry; ++k) {
		__indices_cached(k) = edge.edge(k);
	}
	__vertex_submat.submat(0, 0, last_entry - 1, last_entry - 1) =
		gram.submat(__indices_cached.head(last_entry),
				__indices_cached.head(last_entry));
	for(arma::uword i = 0, max = gram.n_cols; i < max; ++i) {
		if(i == edge.vertex_ind) continue;
		if(std::find(edge.edge.begin(), edge.edge.end(), i) == edge.edge.end()) {
			arma::uvec extra_index = {i};
			__indices_cached(last_entry) = i;
			__vertex_submat.col(last_entry) = gram.submat(__indices_cached, extra_index);
			__vertex_submat.row(last_entry) = gram.submat(extra_index, __indices_cached);
			if(is_elliptic(__vertex_submat)) {
				return i;
			}
		}
	}
	return no_vertex;
}
/* The initial vertex could be any elliptic subdiagram, however the way
 * PolytopeCandidates are constructed mean that the initial gram matrix would be
 * a good choice. The vectors corresponding to this matrix all have zeros in the
 * final coordinate, while any other vectors cannot. Hence this comes down to a
 * search for the vectors which have this zero. */
arma::uvec
PolytopeCheck::initial_vertex(const PolytopeCandidate & p) const {
	const VectorFamily & vf = p.vector_family();
	const arma::uword & last_val = p.real_dimension();
	arma::uvec result(last_val);
	arma::uword result_ind = 0;
	for(arma::uword i = 0, max = vf.size();
			result_ind < last_val && i < max;
			++i) {
		auto vector = vf.unsafe_get(i);
		if(vector(last_val) == 0) {
			result(result_ind) = i;
			++result_ind;
		}
	}
	return result;
}
arma::uvec
PolytopeCheck::construct_edge(const arma::uvec & vertex,
		const arma::uword & ind) const {
	arma::uword length = vertex.size() - 1;
	arma::uvec result(length);
	for(arma::uword i = 0, j = 0; i < length; ++i, ++j) {
		if(j == ind) {
			--i;
			continue;
		}
		result(i) = vertex(j);
	}
	return result;
}
void
PolytopeCheck::add_edges_from_vertex(const arma::uvec & vertex,
		const arma::uword & exclude) {
	for(arma::uword i = 0, max = vertex.size(); i < max; ++i) {
		if(vertex(i) == exclude) continue;
		_edge_queue.emplace(vertex(i), construct_edge(vertex, i));
	}
}
arma::uvec
PolytopeCheck::get_vertex_from_edge(const Edge & cur_edge,
		const arma::uword & vertex_index) const {
	arma::uvec new_vert(cur_edge.edge.size() + 1);
	const arma::uword & last_entry = cur_edge.edge.size();
	for(arma::uword i = 0; i < last_entry; ++i) {
		new_vert(i) = cur_edge.edge(i);
	}
	new_vert(last_entry) = vertex_index;
	std::sort(new_vert.begin(), new_vert.end());
	return new_vert;
}
bool
PolytopeCheck::vertex_unvisited(const arma::uvec & vertex) const {
	return _visited_vertices.find(vertex) == _visited_vertices.end();
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
bool
PolytopeCheck::is_elliptic(const arma::mat & mat) const {
	bool success = has_chol(mat);
	return success;
}
/* Simplifies what arma::chol computes. The upper triangle of
 * __elliptic_check_tmp will contain the actual cholesky decomp of mat. The
 * lower triangle contains junk that would otherwise be set to 0. We don't
 * actually care about the cholesky result here, just whether it can be
 * computed. */
bool
PolytopeCheck::has_chol(const arma::mat & mat) const {
	__elliptic_check_tmp = mat;
	char uplo = 'U';
	arma::blas_int n = __elliptic_check_tmp.n_rows;
	arma::blas_int info = 0;
	arma::lapack::potrf(&uplo, &n, __elliptic_check_tmp.memptr(), &n, &info);
	return (info == 0);
}
}

