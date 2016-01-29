/*
 * polytope_extender_test.cc
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
#include "polytope_extender.h"

#include <gtest/gtest.h>

#include "angles.h"
#include "calc.h"
#include "elliptic_factory.h"

namespace ptope {
TEST(PolytopeExtender, RealCount2) {
	PolytopeCandidate p({{ 1, -.5 }, {-.5, 1} });
	Angles::get().set_angles({2,3});
	PolytopeExtender ext(p);
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(0, count);
	EXPECT_FALSE(ext.has_next());
}
TEST(PolytopeExtender, RealCount3) {
	PolytopeCandidate p({{ 1, -.5, 0 }, {-.5, 1, -.5 }, { 0, -.5, 1 } });
	Angles::get().set_angles({2,3});
	PolytopeExtender ext(p);
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(3, count);
	EXPECT_FALSE(ext.has_next());
}
TEST(PolytopeExtender, HypCount3) {
	PolytopeCandidate p({{ 1, -.5, 0 }, {-.5, 1, -.5 }, { 0, -.5, 1 } });
	PolytopeCandidate q = p.extend_by_inner_products({ -.5, -.5, -.5});
	Angles::get().set_angles({2,3});
	PolytopeExtender ext(q);
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(3, count);
	EXPECT_FALSE(ext.has_next());
}
TEST(PolytopeExtender, NoStop) {
	PolytopeCandidate p(elliptic_factory::type_a(3));
	Angles::get().set_angles({2, 3, 4, 5, 8});
	PolytopeExtender ext(p);
	int count = 0;
	while(ext.has_next()) {
		++count;
		ext.next();
	}
	EXPECT_EQ(115, count);
	EXPECT_FALSE(ext.has_next());
	EXPECT_FALSE(ext.has_next());
	EXPECT_FALSE(ext.has_next());
	EXPECT_FALSE(ext.has_next());
}
/* Stacked Iterator assumes that the first call to next will always work, and
 * does not check has_next. This behaviour is mimicked here. */
TEST(PolytopeExtender, StackedIterMock) {
	using ptope::calc::min_cos_angle;
	Angles::get().set_angles({2, 3, 4, 5, 8, 10});
	PolytopeCandidate p(elliptic_factory::type_a(4));
	auto q = p.extend_by_inner_products({ 0, -.5, min_cos_angle(10), min_cos_angle(10) });
	auto r = q.extend_by_inner_products({ min_cos_angle(10), min_cos_angle(8), 0, 0});
	r.rebase_vectors({ 0, 1, 2, 5 });
	PolytopeExtender ext(r);
	ASSERT_TRUE(ext.next().valid());
}
}

