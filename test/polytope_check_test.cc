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
#include "polytope_check.h"

#include <gtest/gtest.h>

namespace ptope {
namespace {
double min_cos_angle(uint mult) {
	return -std::cos(arma::datum::pi/mult);
}
}
TEST(PolytopeCheck, EsselmannExample) {
	PolytopeCandidate p({ { 1, -.5, 0, 0 }, 
												{ -.5, 1, min_cos_angle(4), 0 }, 
												{ 0, min_cos_angle(4), 1, -.5 }, 
												{ 0, 0, -.5, 1 } });
	auto q = p.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	auto r = q.extend_by_inner_products({ min_cos_angle(8), 0, 0, 0 });
	PolytopeCheck chk;
	EXPECT_TRUE(chk(r));
}
TEST(PolytopeCheck, LannerExample) {
	PolytopeCandidate p({ { 1, -.5, 0 }, { -.5, 1, min_cos_angle(5) },
			{ 0, min_cos_angle(5), 1 }});
	auto q = p.extend_by_inner_products({ 0, 0, -.5});
	PolytopeCheck chk;
	EXPECT_TRUE(chk(q));
}
TEST(PolytopeCheck, Counter) {
	PolytopeCandidate p({ { 1, -.5, 0 }, { -.5, 1, min_cos_angle(5) },
			{ 0, min_cos_angle(5), 1 }});
	auto q = p.extend_by_inner_products({ -.5, -.5, -.5});
	PolytopeCheck chk;
	EXPECT_FALSE(chk(q));
}
}

