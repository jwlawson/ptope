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
bool
PolytopeCheck::operator()(const PolytopeCandidate & p) {
	_visited_vertices.clear();
	while(!_edge_queue.empty()) _edge_queue.pop();

	arma::uvec init_v = initial_vertex(p);
	_visited_vertices.push_back(init_v);
	for(arma::uword i = 0, max = init_v.size(); i < max; ++i) {
		_edge_queue.emplace(i, construct_edge(init_v, i));
	}
	while(!_edge_queue.empty()) {
		Edge cur_edge = _edge_queue.front();
		_edge_queue.pop();
		arma::uword next_vert_ind = find_edge_end(cur_edge, p);
		if(next_vert_ind == no_vertex) {
			return false;
		}
		arma::uvec new_vert = get_vertex_from_edge(cur_edge, next_vert_ind);
		if(vertex_unvisited(new_vert)) {
			_visited_vertices.push_back(new_vert);
			add_edges_from_vertex(new_vert, next_vert_ind);
		}
	}
	return true;
}
arma::uword
PolytopeCheck::find_edge_end(const Edge & edge,
		const PolytopeCandidate & p) const {
	const arma::mat & gram = p.gram();
	arma::uword last_entry = edge.edge.size();
	arma::uword size = last_entry + 1;
	arma::uvec indices(size);
	for(arma::uword k = 0, max = last_entry; k < max; ++k) {
		indices(k) = edge.edge(k);
	}
	for(arma::uword i = 0, max = gram.n_cols; i < max; ++i) {
		if(i == edge.vertex_ind) continue;
		if(std::find(edge.edge.begin(), edge.edge.end(), i) == edge.edge.end()) {
			indices(last_entry) = i;
			auto sub = gram.submat(indices, indices);
			if(is_elliptic(sub)) {
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
	arma::uword last_val = vf.dimension() - 1;
	arma::uvec result(last_val);
	arma::uword result_ind = 0;
	for(arma::uword i = 0, max = vf.size();
			result_ind < last_val && i < max;
			++i) {
		arma::vec vector = vf.unsafe_get(i);
		if(vector(last_val) == 0) {
			result(result_ind) = i;
			++result_ind;
		}
	}
	return result;
}
arma::uvec
PolytopeCheck::construct_edge(const arma::uvec & vertex,
		const arma::uword ind) const {
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
		const arma::uword exclude) {
	_visited_vertices.push_back(vertex);
	for(arma::uword i = 0, max = vertex.size(); i < max; ++i) {
		if(vertex(i) == exclude) continue;
		_edge_queue.emplace(i, construct_edge(vertex, i));
	}
}
arma::uvec
PolytopeCheck::get_vertex_from_edge(const Edge & cur_edge,
		const arma::uword  vertex_index) const {
	arma::uvec new_vert(cur_edge.edge.size() + 1);
	arma::uword last_entry = cur_edge.edge.size();
	for(arma::uword i = 0; i < last_entry; ++i) {
		new_vert(i) = cur_edge.edge(i);
	}
	new_vert(last_entry) = vertex_index;
	std::sort(new_vert.begin(), new_vert.end());
	return new_vert;
}
/* Uses lambda to check whether two uvecs are equal. Cannot use the 'normal'
 * operator== way as that returns a vector rather than a bool. */
bool
PolytopeCheck::vertex_unvisited(const arma::uvec & vertex) const {
	return std::find_if(_visited_vertices.begin(), _visited_vertices.end(),
			[vertex](arma::uvec u) -> bool {arma::uvec diff = vertex - u;
			for(arma::uword val : diff){if(val != 0) return false;} return true;})
		== _visited_vertices.end();
}
bool
PolytopeCheck::is_elliptic(const arma::mat & mat) const {
	if(arma::det(mat) < 0) return false;
	arma::vec evalues(mat.n_cols);
	arma::eig_sym(evalues, mat);
	for(const double & val : evalues) {
		if(val < 0) return false;
	}
	return true;
}
}

