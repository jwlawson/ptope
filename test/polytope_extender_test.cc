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

namespace ptope {
TEST(PolytopeExtender, RealCount2) {
	PolytopeCandidate p({{ 1, -.5 }, {-.5, 1} });
	PolytopeExtender ext(p, {2,3});
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(4, count);
}
TEST(PolytopeExtender, RealCount3) {
	PolytopeCandidate p({{ 1, -.5, 0 }, {-.5, 1, -.5 }, { 0, -.5, 1 } });
	PolytopeExtender ext(p, {2,3});
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(8, count);
}
TEST(PolytopeExtender, HypCount3) {
	PolytopeCandidate p({{ 1, -.5, 0 }, {-.5, 1, -.5 }, { 0, -.5, 1 } });
	PolytopeCandidate q = p.extend_by_inner_products({ -.5, -.5, -.5});
	PolytopeExtender ext(q, {2,3});
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(8, count);
}
}

