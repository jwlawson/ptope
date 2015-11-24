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
#include "elliptic_factory.h"

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
TEST(PolytopeCandidate, StrangeNaN) {
	arma::mat a6 = elliptic_factory::type_a(6);
	PolytopeCandidate p(a6);
	PolytopeCandidate q = p.extend_by_inner_products({ min_cos_angle(5), 0, 0, 0, 0, 0 } );
	std::cout << q.signature().first << " " << q.signature().second << std::endl;
	PolytopeCandidate r = q.extend_by_inner_products({ 0, 0, 0, 0, min_cos_angle(4), min_cos_angle(3) });
	std::cout << r << std::endl;
}
TEST(PolytopeCandidate, SaveLoad) {
	PolytopeCandidate p({ { 1, 0 }, { 0, 1 } });
	p.extend_by_inner_products({ -.5, -.5 });
	std::stringstream ss;
	p.save(ss);
	
	PolytopeCandidate q;
	q.load(ss);

	const arma::mat & p_g = p.gram();
	const arma::mat & q_g = q.gram();
	for(arma::uword i = 0; i < p_g.size(); ++i) {
		EXPECT_DOUBLE_EQ(p_g(i), q_g(i));
	}
}
}

