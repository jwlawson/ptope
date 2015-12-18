/*
 * unique_matrix_check.cc
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
#include "unique_matrix_check.h"

namespace ptope {
bool
UniqueMatrixCheck::operator()(const PolytopeCandidate & p) {
	return operator()(p.gram_ptr());
}
bool
UniqueMatrixCheck::operator()(const std::shared_ptr<const arma::mat> & m) {
	bool result = (_set.find(m) == _set.end());
	if(result) {
		_set.insert(m);
	}
	return result;
}
}
