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
	_initial.rebase_vectors(_perm.next());
	return _initial;
}
PolytopeRebaser::PermIter::PermIter(int size, int max)
	: _next(size),
		_progress(size),
		_progress_max(max),
		_diff(max - size),
		_has_next(true) {
	std::iota(_progress.begin(), _progress.end(), 0);
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
	const std::size_t last = _progress.size() - 1;
	increment(last);
}
void
PolytopeRebaser::PermIter::increment(const std::size_t & ind) {
	if(++_progress[ind] > _diff + ind) {
		if(ind == 0) {
			_has_next = false;
		} else {
			increment(ind - 1);
			_progress[ind] = _progress[ind - 1] + 1;
		}
	}
}
void
PolytopeRebaser::PermIter::compute_next() {
	for(std::size_t i = 0, max = _progress.size(); i < max; ++i) {
		_next(i) = _progress[i];
	}
}
}

