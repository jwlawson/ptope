/*
 * duplicate_column_check_test.cc
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
#include "duplicate_column_check.h"

#include <gtest/gtest.h>

#include "elliptic_factory.h"

namespace ptope {
TEST(DuplicateColumnCheck, Easy) {
	arma::mat a = { { 0, 0 }, { 1, 1 } };
	DuplicateColumnCheck chk;
	EXPECT_TRUE(chk(a));
	arma::mat b = { { 0, 1 }, { 1, 0 } };
	EXPECT_FALSE(chk(b));
}
TEST(DuplicateColumnCheck, Bigger) {
	arma::mat a = { 
	{   1.0000,  -0.8090,  -0.9239,  -0.9239 },
	{  -0.8090,   1.0000,  -0.9239,  -0.9239 },
	{  -0.9239,  -0.9239,   1.0000,   1.0000 },
  {  -0.9239,  -0.9239,   1.0000,   1.0000 } };
	DuplicateColumnCheck chk;
	EXPECT_TRUE(chk(a));
}
}

