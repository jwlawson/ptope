/*
 * compatibility_info_test.cc
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
#include "compatibility_info.h"

#include <gtest/gtest.h>

#include "angles.h"

namespace ptope {
TEST(CompatibilityInfo, NotCompatibleSelf) {
	VectorSet set( 3 );
	set.add( { 1.0, 2.0, 1.0 } );
	set.add( { -4, 0.5, 0.23 } );

	CompatibilityInfo info;
	info.from( set );

	EXPECT_FALSE( info.are_compatible( 0, 0 ) );
	EXPECT_FALSE( info.are_compatible( 1, 1 ) );
}
TEST(CompatibilityInfo, SymmetricFourVectors) {
	VectorSet set( 3, 4 );
	set.add( { 1, 0, 0 } );
	set.add( { 0, 1, 0 } );
	// Note: Should really be using unit vectors, but these provide 'good' inner
	// products, even though they don't correspond to the actual angles.
	set.add( { -0.5, -0.5, 1/std::sqrt(2) } );
	set.add( { std::sqrt(2), 1, std::sqrt(2) } );

	Angles::get().set_angles( { 2, 3 } );
	CompatibilityInfo info;
	info.from( set );

	EXPECT_TRUE( info.are_compatible( 0, 1 ) );
	EXPECT_TRUE( info.are_compatible( 1, 0 ) );
	EXPECT_TRUE( info.are_compatible( 0, 2 ) );
	EXPECT_TRUE( info.are_compatible( 2, 0 ) );
	EXPECT_TRUE( info.are_compatible( 1, 2 ) );
	EXPECT_TRUE( info.are_compatible( 2, 1 ) );
	EXPECT_FALSE( info.are_compatible( 0, 3 ) );
	EXPECT_FALSE( info.are_compatible( 3, 0 ) );
	EXPECT_FALSE( info.are_compatible( 1, 3 ) );
	EXPECT_FALSE( info.are_compatible( 3, 1 ) );
	// <2, 3> gives an inner product smaller than -1
	EXPECT_TRUE( info.are_compatible( 2, 3 ) );
	EXPECT_TRUE( info.are_compatible( 3, 2 ) );
}
TEST(CompatibilityInfo, SymmetricFiveVectors) {
	VectorSet set( 3 );
	set.add( { 1, 0, 0 } );
	set.add( { 0, 1, 0 } );
	set.add( { -1/std::sqrt(2), 1, 1/std::sqrt(2) } );
	set.add( { -0.5, 1, 0.5 } );
	set.add( { 0, 1.5, 1/std::sqrt(2) } );

	Angles::get().set_angles( { 2, 3 } );
	CompatibilityInfo info;
	info.from( set );

	EXPECT_TRUE( info.are_compatible( 0, 1 ) );
	EXPECT_TRUE( info.are_compatible( 1, 0 ) );
	EXPECT_FALSE( info.are_compatible( 0, 2 ) );
	EXPECT_FALSE( info.are_compatible( 2, 0 ) );
	EXPECT_TRUE( info.are_compatible( 0, 3 ) );
	EXPECT_TRUE( info.are_compatible( 3, 0 ) );
	EXPECT_TRUE( info.are_compatible( 0, 4 ) );
	EXPECT_TRUE( info.are_compatible( 4, 0 ) );
	EXPECT_FALSE( info.are_compatible( 1, 2 ) );
	EXPECT_FALSE( info.are_compatible( 2, 1 ) );
	EXPECT_FALSE( info.are_compatible( 1, 3 ) );
	EXPECT_FALSE( info.are_compatible( 3, 1 ) );
	EXPECT_FALSE( info.are_compatible( 1, 4 ) );
	EXPECT_FALSE( info.are_compatible( 4, 1 ) );
	EXPECT_FALSE( info.are_compatible( 2, 3 ) );
	EXPECT_FALSE( info.are_compatible( 3, 2 ) );
	EXPECT_FALSE( info.are_compatible( 2, 4 ) );
	EXPECT_FALSE( info.are_compatible( 4, 2 ) );
}
TEST(CompatibilityInfo, FindNextFiveVectors) {
	VectorSet set( 3 );
	set.add( { 1, 0, 0 } );
	set.add( { 0, 1, 0 } );
	set.add( { std::sqrt(2), 1, std::sqrt(2) } );
	set.add( { -0.5, 1, 0.5 } );
	set.add( { -0.5, 2, std::sqrt(13)/2 } );

	Angles::get().set_angles( { 2, 3 } );
	CompatibilityInfo info;
	info.from( set );

	EXPECT_TRUE( info.are_compatible( 0, 1 ) );
	EXPECT_FALSE( info.are_compatible( 0, 2 ) );
	EXPECT_TRUE( info.are_compatible( 0, 3 ) );
	EXPECT_TRUE( info.are_compatible( 0, 4 ) );

	std::size_t next;
	next = info.next_compatible_to( 0, 0 );
	EXPECT_EQ( 1u, next );
	next = info.next_compatible_to( 0, 1 );
	EXPECT_EQ( 3u, next );
	next = info.next_compatible_to( 0, 3 );
	EXPECT_EQ( 4u, next );

	next = info.next_compatible_to( 0, 4 );
	EXPECT_EQ( 0u, next );

	EXPECT_FALSE( info.are_compatible( 1, 2 ) );
	EXPECT_FALSE( info.are_compatible( 1, 3 ) );
	EXPECT_FALSE( info.are_compatible( 1, 4 ) );
	next = info.next_compatible_to( 1, 0 );
	EXPECT_EQ( 1u, next );

	EXPECT_FALSE( info.are_compatible( 2, 3 ) );
	EXPECT_TRUE( info.are_compatible( 2, 4 ) );
	next = info.next_compatible_to( 2, 0 );
	EXPECT_EQ( 4u, next );
	next = info.next_compatible_to( 2, 3 );
	EXPECT_EQ( 4u, next );
	next = info.next_compatible_to( 2, 4 );
	EXPECT_EQ( 2u, next );
}
TEST(CompatibilityInfo, FindNextSixVectors) {
	VectorSet set( 3 );
	set.add( { 1, 0, 0 } );
	set.add( { 0, 1, 0 } );
	set.add( { std::sqrt(2), 1, std::sqrt(2) } );
	set.add( { -0.5, 1, 0.5 } );
	set.add( { -0.5, 2, std::sqrt(13)/2 } );
	set.add( { 0, 1.5, 1/std::sqrt(2) } );

	Angles::get().set_angles( { 2, 3 } );
	CompatibilityInfo info;
	info.from( set );

	EXPECT_TRUE( info.are_compatible( 0, 1 ) );
	EXPECT_FALSE( info.are_compatible( 0, 2 ) );
	EXPECT_TRUE( info.are_compatible( 0, 3 ) );
	EXPECT_TRUE( info.are_compatible( 0, 4 ) );
	EXPECT_TRUE( info.are_compatible( 0, 5 ) );

	std::size_t next;
	next = info.next_compatible_to( 0, 0 );
	EXPECT_EQ( 1u, next );
	next = info.next_compatible_to( 0, 1 );
	EXPECT_EQ( 3u, next );
	next = info.next_compatible_to( 0, 3 );
	EXPECT_EQ( 4u, next );
	next = info.next_compatible_to( 0, 4 );
	EXPECT_EQ( 5u, next );

	next = info.next_compatible_to( 0, 5 );
	EXPECT_EQ( 0u, next );

	next = info.next_compatible_to( 1, 0 );
	EXPECT_EQ( 1u, next );

	EXPECT_FALSE( info.are_compatible( 2, 3 ) );
	EXPECT_TRUE( info.are_compatible( 2, 4 ) );
	next = info.next_compatible_to( 2, 0 );
	EXPECT_EQ( 4u, next );
	next = info.next_compatible_to( 2, 3 );
	EXPECT_EQ( 4u, next );
	next = info.next_compatible_to( 2, 4 );
	EXPECT_EQ( 2u, next );
}
}

