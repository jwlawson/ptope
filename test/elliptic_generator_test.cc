/*
 * elliptic_generator_test.cc
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
#include "elliptic_generator.h"

#include <gtest/gtest.h>

namespace ptope {
TEST(EllipticGenerator, NumberInSize2) {
	EllipticGenerator gen(2);
	int count = 0;
	while(gen.has_next()) {
		arma::mat a = gen.next();
		EXPECT_EQ((uint) 2, a.n_rows);
		EXPECT_EQ((uint) 2, a.n_cols);
		count++;
	}
	EXPECT_EQ(5, count);
}
TEST(EllipticGenerator, NumberInSize3) {
	EllipticGenerator gen(3);
	int count = 0;
	while(gen.has_next()) {
		arma::mat a = gen.next();
		EXPECT_EQ((uint) 3, a.n_rows);
		EXPECT_EQ((uint) 3, a.n_cols);
		count++;
	}
	EXPECT_EQ(3, count);
}
TEST(EllipticGenerator, NumberInSize6) {
	EllipticGenerator gen(6);
	int count = 0;
	while(gen.has_next()) {
		arma::mat a = gen.next();
		EXPECT_EQ((uint) 6, a.n_rows);
		EXPECT_EQ((uint) 6, a.n_cols);
		count++;
	}
	EXPECT_EQ(4, count);
}
}

