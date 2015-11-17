/*
 * elliptic_generator.cc
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
#include "elliptic_generator.h"

namespace ptope {
namespace {
const std::vector<uint> default_g = { 5, 8, 10};
}
EllipticGenerator::EllipticGenerator(uint size)
	: _size(size) {
	initialize();
}
bool EllipticGenerator::has_next() const {
	return !(todo_g.empty() && done_a && done_b && done_d && done_e && done_f
		&& done_h);
}
arma::mat EllipticGenerator::next() {
	if(!done_a) {
		done_a = true;
		return elliptic_factory::type_a(_size);
	} else if (!done_b) {
		done_b = true;
		return elliptic_factory::type_b(_size);
	} else if (!done_d) {
		done_d = true;
		return elliptic_factory::type_d(_size);
	} else if (!done_e) {
		done_e = true;
		return elliptic_factory::type_e(_size);
	} else if (!done_f) {
		done_f = true;
		return elliptic_factory::type_f(_size);
	} else if (!done_h) {
		done_h = true;
		return elliptic_factory::type_h(_size);
	} else if (!todo_g.empty()) {
		uint label = todo_g.back();
		todo_g.pop_back();
		return elliptic_factory::type_g(_size, label);
	}
	return arma::mat();
}
void EllipticGenerator::initialize() {
	switch(_size) {
		case 2:
			done_d = true;
			todo_g = default_g;
			break;
		case 3:
			done_d = true;
			done_h = false;
			break;
		case 4:
			done_f = false;
			done_h = false;
			break;
		case 6:
		case 7:
		case 8:
			done_e = false;
			break;
		default:
			break;
	}
}
}

