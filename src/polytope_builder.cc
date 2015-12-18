/*
 * polytope_builder.cc
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
#include "polytope_builder.h"

#include "angles.h"

namespace ptope {
PolytopeBuilder::iterator
PolytopeBuilder::build(const PolytopeCandidate & p){
	compute_polytopes(p);
	return _results.begin();
}
PolytopeBuilder::iterator
PolytopeBuilder::end() {
	return _results.end();
}
void
PolytopeBuilder::compute_polytopes(const PolytopeCandidate & p) {
	_results.clear();
	if(_chk(p)) return;
	extend_till_polytope(p, 0);
}
void
PolytopeBuilder::extend_till_polytope(const PolytopeCandidate & p, int depth) {
	if(depth == max_depth) {
		return;
	}
	_vectors[depth].reset();
	auto last_edge = _chk.last_edge();
	auto edge_inv = last_edge.i();
	while(_vectors[depth].has_next()) {
		auto subvec = _ext_vec.subvec(0, _ext_vec.size());
		subvec = _vectors[depth].next();
		if(arma::dot(subvec, edge_inv * subvec) < 1) {
			// elliptic vector!
			for(const double & val : Angles::get().inner_products()) {
				//TODO Check that _ext_vec is actually set properly
				_ext_vec(_ext_vec.size() - 1) = val;
				PolytopeCandidate new_p;
				//Need to rebase the polytope, so that the edge vectors are first!
				//TODO
				bool success = p.extend_by_inner_products(new_p, _ext_vec);
				if(success) {
					//TODO Not sure how to ensure that _chk is in right place to resume
					if(_chk.resume(new_p)){
						_results.push_back(new_p);
					} else {
						extend_till_polytope(new_p, depth + 1);
					}
				}
			}
		}
	}
}
}

