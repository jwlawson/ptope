/*
 * parabolic_filter_test.cc
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
#include "block_diagonal_check.h"

#include <gtest/gtest.h>

namespace ptope {
TEST(BlockDiagonalCheck, Identity2) {
	arma::mat id(2,2);
	id.eye();
	BlockDiagonalCheck filter;
	EXPECT_FALSE(filter(id));
}
TEST(BlockDiagonalCheck, Full2) {
	arma::mat a = { { 1, 2 }, { 3, 4 } };
	BlockDiagonalCheck filter;
	EXPECT_TRUE(filter(a));
}
TEST(BlockDiagonalCheck, LargeDisconnected) {
	arma::mat b = 
	{ { 1, 1, 1, 0, 0, 0, 0 }, 
		{ 1, 1, 1, 0, 0, 0, 0 },
		{ 1, 1, 1, 0, 0, 0, 0 },
		{ 0, 0, 0, 1, 0, 0, 0 },
		{ 0, 0, 0, 0, 1, 1, 0 },
		{ 0, 0, 0, 0, 1, 1, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 } };
	BlockDiagonalCheck filter;
	EXPECT_FALSE(filter(b));
}
}

