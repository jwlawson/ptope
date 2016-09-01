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

	auto first = gram.at_rfp( 0 );
	EXPECT_EQ( 1u, first.row );
	EXPECT_EQ( 1u, first.col );
	EXPECT_DOUBLE_EQ( 2, first.val );

	auto second = gram.at_rfp( 1 );
	EXPECT_EQ( 0u, second.row );
	EXPECT_EQ( 0u, second.col );
	EXPECT_DOUBLE_EQ( 1, second.val );
	
	auto third = gram.at_rfp( 2 );
	EXPECT_EQ( 1u, third.row );
	EXPECT_EQ( 0u, third.col );
	EXPECT_DOUBLE_EQ( 1, third.val );
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
TEST(GramMatrix, Four3dVectors) {
	VectorSet set( 3, 4 );
	set.add( { 2, 0, 0 } );
	set.add( { 0, 3, 0 } );
	set.add( { 1, 1, 1 } );
	set.add( { 2, 1, 2 } );

	arma::vec a = set.at( 0 );
	EXPECT_DOUBLE_EQ( a[0], 2 );
	EXPECT_DOUBLE_EQ( a[1], 0 );
	EXPECT_DOUBLE_EQ( a[2], 0 );
	arma::vec b = set.at( 1 );
	EXPECT_DOUBLE_EQ( b[0], 0 );
	EXPECT_DOUBLE_EQ( b[1], 3 );
	EXPECT_DOUBLE_EQ( b[2], 0 );

	GramMatrix gram;
	gram.from( set );

	EXPECT_DOUBLE_EQ( 4, gram.at( 0, 0 ) );
	EXPECT_DOUBLE_EQ( 9, gram.at( 1, 1 ) );
	EXPECT_DOUBLE_EQ( 1, gram.at( 2, 2 ) );
	EXPECT_DOUBLE_EQ( 1, gram.at( 3, 3 ) );

	EXPECT_DOUBLE_EQ( 0, gram.at( 0, 1 ) );
	EXPECT_DOUBLE_EQ( 0, gram.at( 1, 0 ) );
	EXPECT_DOUBLE_EQ( 2, gram.at( 0, 2 ) );
	EXPECT_DOUBLE_EQ( 2, gram.at( 2, 0 ) );
	EXPECT_DOUBLE_EQ( 4, gram.at( 0, 3 ) );
	EXPECT_DOUBLE_EQ( 4, gram.at( 3, 0 ) );

	EXPECT_DOUBLE_EQ( 3, gram.at( 1, 2 ) );
	EXPECT_DOUBLE_EQ( 3, gram.at( 2, 1 ) );
	EXPECT_DOUBLE_EQ( 3, gram.at( 1, 3 ) );
	EXPECT_DOUBLE_EQ( 3, gram.at( 3, 1 ) );

	EXPECT_DOUBLE_EQ( 1, gram.at( 2, 3 ) );
	EXPECT_DOUBLE_EQ( 1, gram.at( 3, 2 ) );
}
TEST(GramMatrix, Three3dVectors) {
	VectorSet set( 3, 3 );
	set.add( { 2, 0, 1 } );
	set.add( { 0, 3, 0 } );
	set.add( { 2, 1, 2 } );
	GramMatrix gram;
	gram.from( set );

	EXPECT_DOUBLE_EQ( 3, gram.at( 0, 0 ) );
	EXPECT_DOUBLE_EQ( 9, gram.at( 1, 1 ) );
	EXPECT_DOUBLE_EQ( 1, gram.at( 2, 2 ) );

	EXPECT_DOUBLE_EQ( 0, gram.at( 0, 1 ) );
	EXPECT_DOUBLE_EQ( 0, gram.at( 1, 0 ) );
	EXPECT_DOUBLE_EQ( 2, gram.at( 0, 2 ) );
	EXPECT_DOUBLE_EQ( 2, gram.at( 2, 0 ) );

	EXPECT_DOUBLE_EQ( 3, gram.at( 1, 2 ) );
	EXPECT_DOUBLE_EQ( 3, gram.at( 2, 1 ) );
}
TEST(GramMatrix, CheckDiagonalFiveVectors) {
	VectorSet set(3);
	set.add( { 1, 0, 0 } );
	set.add( { 0, 2, 0 } );
	set.add( { 0, 0, 3 } );
	set.add( { 4, 4, 0 } );
	set.add( { 5, 0, 5 } );

	GramMatrix gram;
	gram.from( set );

	EXPECT_DOUBLE_EQ( 1, gram.at( 0, 0 ) );
	EXPECT_DOUBLE_EQ( 4, gram.at( 1, 1 ) );
	EXPECT_DOUBLE_EQ( -9, gram.at( 2, 2 ) );
	EXPECT_DOUBLE_EQ( 32, gram.at( 3, 3 ) );
	EXPECT_DOUBLE_EQ( 0, gram.at( 4, 4 ) );
}
TEST(GramMatrix, CheckDiagonalSixVectors) {
	VectorSet set(3);
	set.add( { 1, 0, 0 } );
	set.add( { 0, 2, 0 } );
	set.add( { 0, 0, 3 } );
	set.add( { 4, 4, 0 } );
	set.add( { 5, 0, 5 } );
	set.add( { 0, 3, 6 } );

	GramMatrix gram;
	gram.from( set );

	EXPECT_DOUBLE_EQ( 1, gram.at( 0, 0 ) );
	EXPECT_DOUBLE_EQ( 4, gram.at( 1, 1 ) );
	EXPECT_DOUBLE_EQ( -9, gram.at( 2, 2 ) );
	EXPECT_DOUBLE_EQ( 32, gram.at( 3, 3 ) );
	EXPECT_DOUBLE_EQ( 0, gram.at( 4, 4 ) );
	EXPECT_DOUBLE_EQ( -27, gram.at( 5, 5 ) );
}
TEST(GramMatrix, RFPIndexFiveVectors) {
	VectorSet set(3);
	set.add( { 1, 0, 0 } );
	set.add( { 0, 2, 0 } );
	set.add( { 0, 0, 3 } );
	set.add( { 4, 4, 0 } );
	set.add( { 5, 0, 5 } );

	GramMatrix gram;
	gram.from( set );

	auto e = gram.at_rfp( 0 );
	EXPECT_EQ( 0u, e.row );
	EXPECT_EQ( 0u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 0, 0 ), e.val );

	e = gram.at_rfp( 3 );
	EXPECT_EQ( 3u, e.row );
	EXPECT_EQ( 0u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 3, 0 ), e.val );

	e = gram.at_rfp( 5 );
	EXPECT_EQ( 3u, e.row );
	EXPECT_EQ( 3u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 3, 3 ), e.val );

	e = gram.at_rfp( 6 );
	EXPECT_EQ( 1u, e.row );
	EXPECT_EQ( 1u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 1, 1 ), e.val );

	e = gram.at_rfp( 9 );
	EXPECT_EQ( 4u, e.row );
	EXPECT_EQ( 1u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 4, 1 ), e.val );

	e = gram.at_rfp( 10 );
	EXPECT_EQ( 4u, e.row );
	EXPECT_EQ( 3u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 4, 3 ), e.val );

	e = gram.at_rfp( 11 );
	EXPECT_EQ( 4u, e.row );
	EXPECT_EQ( 4u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 4, 4 ), e.val );

	e = gram.at_rfp( 12 );
	EXPECT_EQ( 2u, e.row );
	EXPECT_EQ( 2u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 2, 2 ), e.val );
}
TEST(GramMatrix, RFPIndexSixVectors) {
	VectorSet set(3);
	set.add( { 1, 0, 0 } );
	set.add( { 0, 2, 0 } );
	set.add( { 0, 0, 3 } );
	set.add( { 4, 4, 0 } );
	set.add( { 5, 0, 5 } );
	set.add( { 0, 3, 6 } );

	GramMatrix gram;
	gram.from( set );

	auto e = gram.at_rfp( 0 );
	EXPECT_EQ( 3u, e.row );
	EXPECT_EQ( 3u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 3, 3 ), e.val );

	e = gram.at_rfp( 1 );
	EXPECT_EQ( 0u, e.row );
	EXPECT_EQ( 0u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 0, 0 ), e.val );

	e = gram.at_rfp( 6 );
	EXPECT_EQ( 5u, e.row );
	EXPECT_EQ( 0u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 5, 0 ), e.val );

	e = gram.at_rfp( 7 );
	EXPECT_EQ( 4u, e.row );
	EXPECT_EQ( 3u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 4, 3 ), e.val );

	e = gram.at_rfp( 8 );
	EXPECT_EQ( 4u, e.row );
	EXPECT_EQ( 4u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 4, 4 ), e.val );

	e = gram.at_rfp( 9 );
	EXPECT_EQ( 1u, e.row );
	EXPECT_EQ( 1u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 1, 1 ), e.val );

	e = gram.at_rfp( 13 );
	EXPECT_EQ( 5u, e.row );
	EXPECT_EQ( 1u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 5, 1 ), e.val );

	e = gram.at_rfp( 14 );
	EXPECT_EQ( 5u, e.row );
	EXPECT_EQ( 3u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 5, 3 ), e.val );

	e = gram.at_rfp( 15 );
	EXPECT_EQ( 5u, e.row );
	EXPECT_EQ( 4u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 5, 4 ), e.val );

	e = gram.at_rfp( 16 );
	EXPECT_EQ( 5u, e.row );
	EXPECT_EQ( 5u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 5, 5 ), e.val );

	e = gram.at_rfp( 17 );
	EXPECT_EQ( 2u, e.row );
	EXPECT_EQ( 2u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 2, 2 ), e.val );

	e = gram.at_rfp( 20 );
	EXPECT_EQ( 5u, e.row );
	EXPECT_EQ( 2u, e.col );
	EXPECT_DOUBLE_EQ( gram.at( 5, 2 ), e.val );
}
}

