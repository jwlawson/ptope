/*
 * polytope_extender.cc
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
#include "polytope_extender.h"

namespace ptope {
PolytopeExtender::PolytopeExtender(const PolytopeCandidate & initial_polytope)
	:	_initial(initial_polytope),
		_inner_product_vectors(_initial.real_dimension()),
		_computed_next(false) {
	/* Need this call to ensure that the iterator is initialized properly */
	has_next();
}
PolytopeExtender::PolytopeExtender(PolytopeCandidate && initial_polytope)
	:	_initial(initial_polytope),
		_inner_product_vectors(_initial.real_dimension()),
		_computed_next(false) {
	has_next();
}
/**
 * Check whether a subsequent call to next() will return a valid polytope.
 */
bool
PolytopeExtender::has_next(){
	if(!_computed_next) {
		_computed_next = compute_next();
	}
	return _computed_next;
}
/**
 * Compute and return the next extended polytope.
 */
const PolytopeCandidate &
PolytopeExtender::next(){
	_computed_next = false;
	return _next;
}
bool
PolytopeExtender::compute_next() {
	bool success = false;
	while(!success && _inner_product_vectors.has_next()) {
		success = _initial.extend_by_inner_products(_next,
				_inner_product_vectors.next());
	}
	return success;
}
}

