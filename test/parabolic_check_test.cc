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

namespace ptope {
TEST(ParabolicCheck, Identity2) {
	arma::mat a(2,2);
	a.eye();
	ParabolicCheck chk;
	EXPECT_FALSE(chk(a));
}
TEST(ParabolicCheck, SmallParabolic) {
	arma::mat a = { { 1, 0, 0 }, { 0, 2, 0 }, { 0, 0, 0 } };
	ParabolicCheck chk;
	EXPECT_TRUE(chk(a));
}
TEST(ParabolicCheck, NegEValues) {
	arma::mat a = { { 1, 2, 1 }, { 2, 1, 4 }, { 1, 4, 2 } };
	ParabolicCheck chk;
	EXPECT_FALSE(chk(a));
}
TEST(ParabolicCheck, PositiveDefinite) {
	arma::mat a = { { 4, 1 }, { 1, 2 } };
	ParabolicCheck chk;
	EXPECT_FALSE(chk(a));
}
}

