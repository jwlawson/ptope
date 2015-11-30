/*
 * polytope_rebaser.cc
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
#include "polytope_rebaser.h"

namespace ptope {
PolytopeRebaser::PolytopeRebaser(const PolytopeCandidate & p)
	: _initial(p), 
		_basis_ind(0),
		_basis_max(p.real_dimension()),
		_swap_ind(_basis_max),
		_swap_max(p.gram().n_cols) {}
bool
PolytopeRebaser::has_next() {
	return !(_swap_ind == _swap_max && _basis_ind == (_basis_max - 1) );
}
PolytopeCandidate
PolytopeRebaser::next() {
	if(_swap_ind == _swap_max) {
		++_basis_ind;
		_swap_ind = _basis_max;
	}
	return _initial.swap_rebase(_basis_ind, _swap_ind++);
}
}

