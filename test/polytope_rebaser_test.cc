/*
 * polytope_rebaser_test.cc
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
#include "polytope_rebaser.h"
#include "elliptic_factory.h"

#include <gtest/gtest.h>

namespace ptope {
TEST(PolytopeRebaser, Simple) {
	PolytopeCandidate p(elliptic_factory::type_b(2));
	PolytopeCandidate r = p.extend_by_inner_products({ -.5, -.5 });
	ASSERT_TRUE(r.valid());
	PolytopeRebaser rebaser(r);
	int count = 0;
	while(rebaser.has_next()) {
		rebaser.next();
		++count;
	}
	ASSERT_EQ(3, count);
	PolytopeRebaser rebaser2(r);
	/* First one out is the same as initial r */
	rebaser2.next();
	ASSERT_TRUE(rebaser2.has_next());
	ASSERT_TRUE(rebaser2.has_next());
	arma::mat g1 = rebaser2.next().gram();
	ASSERT_TRUE(rebaser2.has_next());
	arma::mat g2 = rebaser2.next().gram();
	EXPECT_FALSE(rebaser2.has_next());

	arma::mat exp = {
		{ 1.0000, -0.5000, -0.7071 },
  	{-0.5000,  1.0000, -0.5000 },
  	{-0.7071, -0.5000,  1.0000 } };
	arma::mat diff1 = g1 - exp;
	for(const double & val : diff1) {
		EXPECT_NEAR(0, val, 1e-4);
	}
	arma::mat diff2 = g2 - exp;
	for(const double & val : diff2) {
		EXPECT_NEAR(0, val, 1e-4);
	}
}
}

