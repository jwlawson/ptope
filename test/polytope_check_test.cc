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
#include "elliptic_factory.h"

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
TEST(PolytopeCheck, TumarkinExample) {
	PolytopeCandidate p( { { 1, min_cos_angle(4), 0, 0 },
			{ min_cos_angle(4), 1, -.5, 0 },
			{0, -.5, 1, -.5 },
			{ 0, 0, -.5, 1 } });
	auto q = p.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });
	auto r = q.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
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
TEST(PolytopeCheck, Odd4) {
	PolytopeCheck chk;
	PolytopeCandidate p(elliptic_factory::type_a(4));
	auto q = p.extend_by_inner_products({ min_cos_angle(5), -.5, 0, 0});
	ASSERT_TRUE(q.valid());
	std::cout << q;
	EXPECT_FALSE(chk(q));
	auto r = q.extend_by_inner_products({0, -.5, min_cos_angle(5), -.5 });
	ASSERT_TRUE(r.valid());
	std::cout << r;
	EXPECT_FALSE(chk(r));
	auto s = r.swap_rebase(3,5);
	auto t = s.extend_by_inner_products({ 0, -.5, 0, 0 });
	std::cout << t;
	ASSERT_TRUE(t.valid());
	EXPECT_TRUE(chk(t));
	EXPECT_TRUE(chk.used_all());
	auto u = t.swap_rebase(2,4);
	EXPECT_TRUE(chk(u));
	EXPECT_TRUE(chk.used_all());
	auto v = t.swap_rebase(0, 5);
	std::cout << v;
	EXPECT_TRUE(chk(v));
	EXPECT_TRUE(chk.used_all());
}
}

