/*
 * number_dotted_check_test.cc
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
#include "number_dotted_check.h"
#include "elliptic_factory.h"

#include <gtest/gtest.h>

namespace ptope {
TEST(NumberDottedCheck, NoDots) {
	PolytopeCandidate p(elliptic_factory::type_a(4));
	NumberDottedCheck<0> zero;
	NumberDottedCheck<1> one;
	EXPECT_TRUE(zero(p));
	EXPECT_FALSE(one(p));
}
}

