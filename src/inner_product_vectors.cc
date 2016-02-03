/*
 * inner_product_vectors.cc
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
#include "inner_product_vectors.h"

#include "angles.h"

namespace ptope {
InnerProductVectors::InnerProductVectors(int size)
	:	_inner_products(Angles::get().inner_products()),
		_next(size),
		_progress(size),
		_progress_max(_inner_products.size() - 1),
		_has_next(true) {
	std::reverse(_inner_products.begin(), _inner_products.end());
}
bool
InnerProductVectors::has_next() {
	return _has_next;
}
const arma::vec &
InnerProductVectors::next() {
	compute_next();
	increment_progress();
	return _next;
}
void
InnerProductVectors::reset() {
	_progress.assign(_progress.size(), 0);
	_has_next = true;
}
void
InnerProductVectors::increment_progress() {
	bool rollover = true;
	std::size_t i = 0;
	const std::size_t max = _progress.size();
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
void
InnerProductVectors::compute_next() {
	for(std::size_t i = 0, max = _progress.size(); i < max; ++i) {
		_next(i) = _inner_products[_progress[i]];
	}
}
}

