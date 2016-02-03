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
		_perm(p.real_dimension(), p.gram().n_cols) {}
bool
PolytopeRebaser::has_next() {
	return _perm.has_next();
}
const PolytopeCandidate &
PolytopeRebaser::next() {
	/* TODO This writes the memory of _next twice. Should be optimized. */
	_next = _initial;
	auto & p = _perm.next();
	_next.rebase_vectors(p);
	return _next;
}
PolytopeRebaser::PermIter::PermIter(int size, int max)
	: _next(size),
		_progress(size),
		_last_ind(size - 1),
		_has_next(true) {
	std::iota(_progress.begin(), _progress.end(), max - size);
}
bool
PolytopeRebaser::PermIter::has_next() {
	return _has_next;
}
const arma::uvec &
PolytopeRebaser::PermIter::next() {
	compute_next();
	increment_progress();
	return _next;
}
void
PolytopeRebaser::PermIter::increment_progress() {
	increment(0);
}
void
PolytopeRebaser::PermIter::increment(const std::size_t & ind) {
	if(_progress[ind] == ind) {
		if(ind == _last_ind) {
			_has_next = false;
		} else {
			increment(ind + 1);
			_progress[ind] = _progress[ind + 1] - 1;
		}
	} else {
		--_progress[ind];
	}
}
void
PolytopeRebaser::PermIter::compute_next() {
	for(std::size_t i = 0, max = _progress.size(); i < max; ++i) {
		_next(i) = _progress[i];
	}
}
}

