/*
 * vector_family_test.cc
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
#include "vector_family.h"

#include <gtest/gtest.h>

namespace ptope {
TEST(VectorFamily, Construct2) {
	arma::mat a = { { 1, 2 }, { 0, 1 } };
	VectorFamily v(a);
	arma::vec a0 = { 1, 0 };
	arma::vec v0 = v.get(0);
	EXPECT_DOUBLE_EQ(a0(0), v0(0));
	EXPECT_DOUBLE_EQ(a0(1), v0(1));
	arma::vec a1 = { 2, 1 };
	arma::vec v1 = v.get(1);
	EXPECT_DOUBLE_EQ(a1(0), v1(0));
	EXPECT_DOUBLE_EQ(a1(1), v1(1));
}
TEST(VectorFamily, Add) {
	arma::mat a = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	arma::vec b = { 2, 3, 4 };
	VectorFamily v(a);
	v.add_vector(b);
	arma::vec v3 = v.get(3);
	EXPECT_DOUBLE_EQ(b(0), v3(0));
	EXPECT_DOUBLE_EQ(b(1), v3(1));
	EXPECT_DOUBLE_EQ(b(2), v3(2));
}
TEST(VectorFamily, MoveConstruct) {
	VectorFamily v({ { 1, 0 }, { 0, 1 } });
	arma::vec v0 = v.get(0);
	EXPECT_DOUBLE_EQ(1, v0(0));
	EXPECT_DOUBLE_EQ(0, v0(1));
	arma::vec v1 = v.get(1);
	EXPECT_DOUBLE_EQ(0, v1(0));
	EXPECT_DOUBLE_EQ(1, v1(1));
}
TEST(VectorFamily, GetSafe) {
	VectorFamily v({ { 1, 2 }, { 3, 4 } });
	arma::vec v0 = v.get(0);
	v0(0) = 5;
	arma::vec v1 = v.get(0);
	EXPECT_DOUBLE_EQ(1, v1(0));
}
TEST(VectorFamily, UnsafeGetUnsafe) {
	VectorFamily v({ { 1, 2 }, { 3, 4 } });
	arma::vec v0 = v.unsafe_get(0);
	v0(0) = 5;
	arma::vec v1 = v.get(0);
	EXPECT_DOUBLE_EQ(5, v1(0));
}
}
