/*
 * angle_check.cc
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
#include "angle_check.h"

#include <gtest/gtest.h>

#include "angles.h"

namespace ptope {
namespace {
double min_cos_angle(uint mult) {
	return -std::cos(arma::datum::pi/mult);
}
}
TEST(AngleCheck, Simple2) {
	Angles::get().set_angles({2, 3, 4, 5, 8});
	arma::mat a = { { 1, 0 }, { 0, 1 } };
	AngleCheck chk;
	EXPECT_TRUE(chk(a));

	arma::mat b = { { 1, -.9 }, { -.9, 1 } };
	EXPECT_FALSE(chk(b));
}
TEST(AngleCheck, Simple3) {
	Angles::get().set_angles({2, 3, 4});
	arma::mat a = { { 1, min_cos_angle(4), min_cos_angle(3) },
		{ min_cos_angle(4), 1, min_cos_angle(4) },
		{ min_cos_angle(3), min_cos_angle(4), 1} };
	AngleCheck chk;
	EXPECT_TRUE(chk(a));

	/* a(0) = -1 would return true, as only last column checked */
	/* This is like an invalid gram matrix. */
	a(6) = 1.5;
	EXPECT_FALSE(chk(a));
	/* This is like an angle not in the list. */
	a(6) = -0.9;
	EXPECT_FALSE(chk(a));
	/* This is like a dotted edge. */
	a(6) = -1.9;
	EXPECT_TRUE(chk(a));
}
TEST(AngleCheck, CustomAngles) {
	Angles::get().set_angles({2,3});
	AngleCheck chk;
	arma::mat a = { { 1, 0 }, { 0, 1 } };
	EXPECT_TRUE(chk(a));
	a(2) = min_cos_angle(3);
	EXPECT_TRUE(chk(a));
	a(2) = min_cos_angle(4);
	EXPECT_FALSE(chk(a));
}
TEST(AngleCheck, DottedEdges) {
	Angles::get().set_angles({2, 3, 4, 5, 8});
	arma::mat a = { { 1, min_cos_angle(4), -1.5 },
		{ min_cos_angle(3), 1, min_cos_angle(5) },
		{ min_cos_angle(3), min_cos_angle(4), 1} };
	AngleCheck chk;
	EXPECT_TRUE(chk(a));
}
TEST( AngleCheck, ZeroVal ) {
	Angles::get().set_angles({ 2 });
	AngleCheck chk;
	EXPECT_TRUE( chk(0) );
	EXPECT_TRUE( chk(0.0) );
	EXPECT_TRUE( chk(1e-14) );
	EXPECT_TRUE( chk(-1e-14) );
}
TEST( AngleCheck, LessThanMinusOne ) {
	Angles::get().set_angles({ 2 , 3 , 4 });
	AngleCheck chk;
	EXPECT_TRUE( chk(-1.1) );
	EXPECT_TRUE( chk(-1.0 - 1e-14) );
	EXPECT_TRUE( chk(-1.0 + 1e-19) );
	EXPECT_TRUE( chk(-29) );
}

}

