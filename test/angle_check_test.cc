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

namespace ptope {
namespace {
double min_cos_angle(uint mult) {
	return -std::cos(arma::datum::pi/mult);
}
}
TEST(AngleCheck, Simple2) {
	arma::mat a = { { 1, 0 }, { 0, 1 } };
	AngleCheck chk;
	EXPECT_TRUE(chk(a));

	arma::mat b = { { -.9, 0 }, { 0, -.9 } };
	EXPECT_FALSE(chk(b));
}
TEST(AngleCheck, Simple3) {
	arma::mat a = { { min_cos_angle(3), min_cos_angle(4), min_cos_angle(5) },
		{ min_cos_angle(3), min_cos_angle(4), min_cos_angle(5) },
		{ min_cos_angle(3), min_cos_angle(4), min_cos_angle(5) } };
	AngleCheck chk;
	EXPECT_TRUE(chk(a));

	/* a(0) = -1 would return true, as only last column checked */
	a(8) = -1;
	EXPECT_FALSE(chk(a));
}
TEST(AngleCheck, CustomAngles) {
	AngleCheck chk({2, 3});
	arma::mat a = { { 1, 0 }, { 0, 1 } };
	EXPECT_TRUE(chk(a));
	a(2) = min_cos_angle(3);
	EXPECT_TRUE(chk(a));
	a(3) = min_cos_angle(4);
	EXPECT_FALSE(chk(a));
}
TEST(AngleCheck, DottedEdges) {
	arma::mat a = { { min_cos_angle(3), min_cos_angle(4), 1.5 },
		{ min_cos_angle(3), min_cos_angle(4), min_cos_angle(5) },
		{ min_cos_angle(3), min_cos_angle(4), min_cos_angle(5) } };
	AngleCheck chk;
	EXPECT_TRUE(chk(a));
}
}

