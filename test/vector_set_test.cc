/*
 * vector_set_test.cc
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
#include "vector_set.h"
#include "gtest/gtest.h"

namespace ptope {
TEST(VectorSet, AddSame) {
	VectorSet<double> set( 5 );
	arma::vec a { 1.0, 2.0, 3.0, 4.0, 5.0 };

	EXPECT_FALSE( set.contains(a) );
	ASSERT_TRUE( set.add(a) );
	EXPECT_TRUE( set.contains(a) );
	EXPECT_FALSE( set.add(a) );
}
TEST(VectorSet, ForceResize) {
	VectorSet<double> set( 3 , 2 );
	arma::vec a { 1.0, 1.1, 1.2 };
	arma::vec b { 1.1, 1.2, 1.3 };
	arma::vec c { 1.2, 1.3, 1.4 };

	EXPECT_FALSE( set.contains(a) );
	EXPECT_TRUE( set.add(a) );
	EXPECT_FALSE( set.contains(b) );
	EXPECT_TRUE( set.add(b) );
	EXPECT_FALSE( set.contains(c) );
	// Adding a third vector forces a resize
	EXPECT_TRUE( set.add(c) );
	EXPECT_TRUE( set.contains(a) );
	EXPECT_TRUE( set.contains(b) );

	arma::vec vec = set.at( 0 );
	EXPECT_DOUBLE_EQ( 1.0 , vec[0] );
	EXPECT_DOUBLE_EQ( 1.1 , vec[1] );
	EXPECT_DOUBLE_EQ( 1.2 , vec[2] );
	vec = set.at( 1 );
	EXPECT_DOUBLE_EQ( 1.1 , vec[0] );
	EXPECT_DOUBLE_EQ( 1.2 , vec[1] );
	EXPECT_DOUBLE_EQ( 1.3 , vec[2] );
	vec = set.at( 2 );
	EXPECT_DOUBLE_EQ( 1.2 , vec[0] );
	EXPECT_DOUBLE_EQ( 1.3 , vec[1] );
	EXPECT_DOUBLE_EQ( 1.4 , vec[2] );
}
TEST(VectorSet, CreateIterator) {
	VectorSet<double> set( 3 );
	arma::vec a { 1.0, 1.1, 1.2 };
	arma::vec b { 1.1, 1.2, 1.3 };
	arma::vec c { 1.2, 1.3, 1.4 };
	set.add(a);
	set.add(b);
	set.add(c);

	VectorSet<double>::iterator it = set.begin();
	arma::vec& first = *it;
	EXPECT_DOUBLE_EQ( a[0], first[0] );
	EXPECT_DOUBLE_EQ( a[1], first[1] );
	EXPECT_DOUBLE_EQ( a[2], first[2] );

	++it;
	arma::vec& second = *it;
	EXPECT_DOUBLE_EQ( b[0], second[0] );
	EXPECT_DOUBLE_EQ( b[1], second[1] );
	EXPECT_DOUBLE_EQ( b[2], second[2] );
}
}

