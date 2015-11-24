/*
 * polytope_candidate_test.cc
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
#include "polytope_candidate.h"

#include <gtest/gtest.h>

namespace ptope {
namespace {
double min_cos_angle(uint mult) {
	return -std::cos(arma::datum::pi/mult);
}
}
TEST(PolytopeCandidate, Print) {
	arma::mat gram = { { 1, -.5, 0 }, {-.5, 1, -.5 }, { 0, -.5, 1 } };
	PolytopeCandidate p(gram);
	//std::cout << p << std::endl;

	PolytopeCandidate ext = p.extend_by_inner_products({ min_cos_angle(3),min_cos_angle(3), min_cos_angle(4) });
	//std::cout << ext << std::endl;

	PolytopeCandidate e2 = ext.extend_by_inner_products({ min_cos_angle(3),min_cos_angle(5), min_cos_angle(4) });
}
}

