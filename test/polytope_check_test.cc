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

#include "angles.h"
#include "calc.h"
#include "elliptic_factory.h"
#include "parabolic_check.h"

#include <gtest/gtest.h>

namespace ptope {
using ptope::calc::min_cos_angle;
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
TEST(PolytopeCheck, TumarkinExample) {
	PolytopeCandidate p( { { 1, min_cos_angle(4), 0, 0 },
			{ min_cos_angle(4), 1, -.5, 0 },
			{0, -.5, 1, -.5 },
			{ 0, 0, -.5, 1 } });
	auto q = p.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });
	ASSERT_TRUE(q.valid());
	auto r = q.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	ASSERT_TRUE(r.valid());
	r.rebase_vectors({ 1, 2, 3, 4 });
	auto s = r.extend_by_inner_products({ 0, 0, min_cos_angle(4), 0 });
	ASSERT_TRUE(s.valid());
	PolytopeCheck chk;
	EXPECT_TRUE(chk(s));
}
TEST(PolytopeCheck, Biggest) {
	PolytopeCandidate p(elliptic_factory::type_e(8));
	auto q = p.extend_by_inner_products({ min_cos_angle(5), 0, 0, 0, 0, 0, 0, 0 });
	auto r = q.extend_by_inner_products({ 0, 0, 0, 0, 0, 0, 0, min_cos_angle(5) });
	r.rebase_vectors({ 0, 1, 2, 4, 5, 6, 7, 8 });
	auto s = r.extend_by_inner_products({ 0, 0, 0, 0, -.5, 0, 0, 0 });
	ASSERT_TRUE(s.valid());
	PolytopeCheck chk;
	EXPECT_TRUE(chk(s));

	s.rebase_vectors({ 0, 1, 2, 8, 3, 4, 5, 6});
	EXPECT_TRUE(chk(s));
}
TEST(PolytopeCheck, 9dim) {
	PolytopeCandidate p(elliptic_factory::type_a(9));
	PolytopeCheck chk;
	auto q = p.extend_by_inner_products({ min_cos_angle(5), 0, 0, 0, 0, 0, 0, 0, 0 });
	ASSERT_TRUE(q.valid());
	auto r = q.extend_by_inner_products({0, 0, min_cos_angle(5), 0, -.5, 0, min_cos_angle(5), min_cos_angle(5), 0});
	ASSERT_TRUE(r.valid());
	auto s = r.extend_by_inner_products({ 0, 0, 0, -.5, min_cos_angle(5), 0, 0, 0, -.5 });
	ASSERT_TRUE(s.valid());
	EXPECT_FALSE(chk(s));
	auto m = s.swap_rebase(0, 10);
	EXPECT_FALSE(chk(m));
	auto t = r.swap_rebase(0, 10);
	auto u = t.extend_by_inner_products({ -.5, 0, 0, -.5, min_cos_angle(5), 0, 0, 0, -.5 });
	ASSERT_TRUE(u.valid());
	EXPECT_FALSE(chk(u));
}
/*
 * This test looks at a polytope which my algorithm decided was compact when
 * checked with the rows in a certain order, but not in a different order.
 */
TEST(PolytopeCheck, Odd4) {
	PolytopeCheck chk;
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
	auto u = t.swap_rebase(2,4);
	EXPECT_FALSE(chk(u));
	auto v = t.swap_rebase(0, 5);
	EXPECT_FALSE(chk(v));
}
TEST(PolytopeCheck, Unfolded) {
	PolytopeCheck chk;
	PolytopeCandidate p(elliptic_factory::type_b(4));
	auto q = p.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });
	ASSERT_TRUE(q.valid());
	EXPECT_FALSE(chk(q));
	auto r = q.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	ASSERT_TRUE(r.valid());
	EXPECT_FALSE(chk(r));
	auto s = r.swap_rebase(0, 5);
	auto t = s.extend_by_inner_products({ min_cos_angle(8), 0, -.5, 0 });
	ASSERT_TRUE(t.valid());
	EXPECT_FALSE(chk(t));

	auto u = t.swap_rebase(3, 4);
	auto v = u.extend_by_inner_products({ 0, min_cos_angle(4), 0, 0 });
	ASSERT_TRUE(v.valid());
	EXPECT_TRUE(chk(v));

	r.rebase_vectors({ 1, 2, 4, 5 });
	auto x = r.extend_by_inner_products({ min_cos_angle(4), 0, 0, 0 });
	ASSERT_TRUE(x.valid());
	auto y = x.extend_by_inner_products({ 0, -.5, 0, min_cos_angle(8) });
	ASSERT_TRUE(y.valid());
	EXPECT_TRUE(chk(v));
}
TEST(PolytopeCheck, 4Dim) {
	PolytopeCheck chk;
	Angles::get().set_angles({2, 3, 4, 5, 8, 10});
	PolytopeCandidate p(elliptic_factory::type_b(4));
	auto q = p.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });
	ASSERT_TRUE(q.valid());
	ASSERT_FALSE(chk(q));
	auto r = q.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	ASSERT_TRUE(r.valid());
	ASSERT_FALSE(chk(r));

	PolytopeCandidate copy(r);
	copy.rebase_vectors({ 1, 2, 4, 5});
	auto c1 = copy.extend_by_inner_products({ min_cos_angle(4), 0, 0, 0 });
	ASSERT_TRUE(c1.valid());
	ASSERT_FALSE(chk(c1));
	auto c2 = c1.extend_by_inner_products({ 0, -.5, 0, min_cos_angle(8) });
	ASSERT_TRUE(c2.valid());
	ParabolicCheck para;
	EXPECT_FALSE(para(c2));
	EXPECT_TRUE(chk(c2));
}
}

