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
namespace {
std::vector<double> angles_to_prods(const std::vector<uint> & angles) {
	std::vector<double> result(angles.size());
	for(std::size_t i = 0, max = angles.size(); i < max; ++i) {
		if(angles[i] == 2) {
			result[i] = 0;
		} else {
			result[i] = -std::cos(arma::datum::pi/angles[i]);
		}
	}
	return result;
}
}
std::vector<double>
PolytopeExtender::__default_inner = angles_to_prods({ 2, 3, 4, 5, 8 });
PolytopeExtender::PolytopeExtender(const PolytopeCandidate & initial_polytope)
	:	_initial(initial_polytope),
		_inner_products(__default_inner),
		_extend_by(_initial.real_dimension()),
		_progress(_initial.real_dimension()),
		_progress_max(_inner_products.size() - 1),
		_has_next(true) {}
PolytopeExtender::PolytopeExtender(PolytopeCandidate && initial_polytope)
	:	_initial(initial_polytope),
		_inner_products(__default_inner),
		_extend_by(_initial.real_dimension()),
		_progress(_initial.real_dimension()),
		_progress_max(_inner_products.size() - 1),
		_has_next(true) {}
PolytopeExtender::PolytopeExtender(const PolytopeCandidate & initial_polytope,
		const std::vector<uint> & angles)
	:	_initial(initial_polytope),
		_inner_products(angles_to_prods(angles)),
		_extend_by(_initial.real_dimension()),
		_progress(_initial.real_dimension()),
		_progress_max(_inner_products.size() - 1),
		_has_next(true) {}
PolytopeExtender::PolytopeExtender(PolytopeCandidate && initial_polytope,
		const std::vector<uint> & angles)
	:	_initial(initial_polytope),
		_inner_products(angles_to_prods(angles)),
		_extend_by(_initial.real_dimension()),
		_progress(_initial.real_dimension()),
		_progress_max(_inner_products.size() - 1),
		_has_next(true) {}
/**
 * Check whether a subsequent call to next() will return a valid polytope.
 */
bool
PolytopeExtender::has_next(){
	return _has_next;
}
/**
 * Compute and return the next extended polytope.
 */
PolytopeCandidate
PolytopeExtender::next(){
	for(std::size_t i = 0, max = _progress.size(); i < max; ++i) {
		_extend_by(i) = _inner_products[_progress[i]];
	}
	PolytopeCandidate result = _initial.extend_by_inner_products(_extend_by);	
	increment_progress();
	return std::move(result);
}
void
PolytopeExtender::set_default_angles(const std::vector<uint> & angles) {
	__default_inner = angles_to_prods(angles);
}
void
PolytopeExtender::increment_progress() {
	bool rollover = true;
	std::size_t i = 0;
	std::size_t max = _progress.size();
	for(; rollover && i < max; ++i) {
		rollover = false;
		if(++_progress[i] > _progress_max) {
			_progress[i] = 0;
			rollover = true;
		}
	}
	if(i == max && rollover) {
		_has_next = false;
	}
}
}

