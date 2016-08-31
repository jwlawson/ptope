/*
 * gram_matrix.cc
 * Copyright 2016 John Lawson
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
#include "gram_matrix.h"
#include "gtest/gtest.h"

namespace ptope {
TEST(GramMatrix, TwoEuclVectors) {
	VectorSet set(3, 2);
	set.add( { 1, 0, 0 } );
	set.add( { 1, 1, 0 } );
	GramMatrix gram;
	gram.from( set );
	double aa = gram.at( 0, 0 );
	double ab = gram.at( 0, 1 );
	double ba = gram.at( 1, 0 );
	double bb = gram.at( 1, 1 );

	EXPECT_DOUBLE_EQ( 1, aa );
	EXPECT_DOUBLE_EQ( 1, ab );
	EXPECT_DOUBLE_EQ( 1, ba );
	EXPECT_DOUBLE_EQ( 2, bb );

	auto first = gram.at( 0 );
	EXPECT_EQ( 0, first.row );
	EXPECT_EQ( 1, first.col );
	EXPECT_EQ( 1, first.val );

	auto second = gram.at( 1 );
	EXPECT_EQ( 1, second.row );
	EXPECT_EQ( 1, second.col );
	EXPECT_EQ( 2, second.val );
	
	auto third = gram.at( 2 );
	EXPECT_EQ( 0, third.row );
	EXPECT_EQ( 0, third.col );
	EXPECT_EQ( 1, third.val );
}
TEST(GramMatrix, TwoHypVectors) {
	VectorSet set(3, 2);
	set.add( { 1, 0, 1 } );
	set.add( { 0, 1, 2 } );
	GramMatrix gram;
	gram.from( set );
	double aa = gram.at( 0, 0 );
	double ab = gram.at( 0, 1 );
	double ba = gram.at( 1, 0 );
	double bb = gram.at( 1, 1 );

	EXPECT_DOUBLE_EQ( 0, aa );
	EXPECT_DOUBLE_EQ( -2, ab );
	EXPECT_DOUBLE_EQ( -2, ba );
	EXPECT_DOUBLE_EQ( -3, bb );
}
}

