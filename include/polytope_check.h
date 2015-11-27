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

#include "polytope_candidate.h"

#include <queue>

namespace ptope {
class PolytopeCheck {
private:
	static constexpr arma::uword
	no_vertex = std::numeric_limits<arma::uword>::max();
	struct Edge {
		Edge(arma::uword & ind, arma::uvec && vec) 
			: vertex_ind(ind), edge(vec) {}
		arma::uword vertex_ind;
		arma::uvec edge;
	};
public:
	/**
	 * Check whether a given polytope candidate is actually a compact polytope.
	 * Returns true if it is a compact polytope.
	 */
	bool operator()(const PolytopeCandidate & p);
private:
	std::vector<arma::uvec> _visited_vertices;
	std::queue<Edge> _edge_queue;
	/**
	 * Find the vertex at the end of an edge. Each edge is constructed from an
	 * initial vertex, so this finds the other vertex along the edge. If no vertex
	 * is found then PolytopeCheck::no_vertex is returned.
	 *
	 * The index returned is the additional vertex index needed to complete the
	 * edge to a vertex.
	 */
	arma::uword
	find_edge_end(const Edge & edge, const PolytopeCandidate & p) const;
	/**
	 * Find an initial elliptic subdiagram to use as initial vertex.
	 */
	arma::uvec
	initial_vertex(const PolytopeCandidate & p) const;
	/**
	 * Construct an edge adjacent to the given vertex by removing the specified
	 * index.
	 */
	arma::uvec
	construct_edge(const arma::uvec & vertex, const arma::uword index) const;
	/**
	 * Use the provided vertes to add all edges adjacent to that vertex to the
	 * queue, except for the edge specified by the excluded index. The excluded
	 * index corresponds to the edge which first visited this vertex.
	 */
	void
	add_edges_from_vertex(const arma::uvec & vertex,
			const arma::uword exclude);
	/**
	 * Construct the vector representing the vertex which lies at the end of the
	 * specified edge with the given additional index.
	 */
	arma::uvec
	get_vertex_from_edge(const Edge & cur_edge,
			const arma::uword vertex_index) const;
	/**
	 * Check whether the vertex has been visited before.
	 */
	bool
	vertex_unvisited(const arma::uvec & vertex) const;
	/**
	 * Check whether the given matrix is elliptic (i.e. positive definite).
	 */
	bool
	is_elliptic(const arma::mat & mat) const;
};
}
#endif
