/*
 * parabolic_check.h
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
#include "parabolic_check.h"

#include <gtest/gtest.h>

#include "elliptic_factory.h"

namespace ptope {
namespace {
double min_cos_angle(uint mult) {
	return -std::cos(arma::datum::pi/mult);
}
}
TEST(ParabolicCheck, QuasiLanner) {
	PolytopeCandidate p(elliptic_factory::type_b(4));
	ParabolicCheck chk;
	EXPECT_FALSE(chk(p));
	auto q = p.extend_by_inner_products({ -.5, 0, 0, -std::sqrt(2)/2 });
	ASSERT_TRUE(q.valid());
	EXPECT_TRUE(chk(q));
}
/*
 * The following test looks at the same 4+3 dimensional polytope as the
 * PolytopeCheck.Odd4 test. Here we see that at no point does this contain any
 * parabolics.
 */
TEST(ParabolicCheck, 4Plus3) {
	ParabolicCheck chk;
	PolytopeCandidate p(elliptic_factory::type_a(4));
	auto q = p.extend_by_inner_products({ min_cos_angle(5), -.5, 0, 0});
	ASSERT_TRUE(q.valid());
	EXPECT_FALSE(chk(q));
	auto r = q.extend_by_inner_products({0, -.5, min_cos_angle(5), -.5 });
	ASSERT_TRUE(r.valid());
	EXPECT_FALSE(chk(r));
	auto s = r.swap_rebase(3,5);
	auto t = s.extend_by_inner_products({ 0, -.5, 0, 0 });
	ASSERT_TRUE(t.valid());
	EXPECT_FALSE(chk(t));
}
}

