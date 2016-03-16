/*
 * polytope_extender_test.cc
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
#include "polytope_extender.h"

#include <gtest/gtest.h>

#include "angles.h"
#include "calc.h"
#include "elliptic_factory.h"
#include "stacked_iterator.h"
#include "combined_check.h"
#include "polytope_rebaser.h"
#include "angle_check.h"
#include "filtered_iterator.h"
#include "duplicate_column_check.h"
#include "parabolic_check.h"
#include "unique_matrix_check.h"
#include "construct_iterator.h"
#include "polytope_check.h"
#include "elliptic_generator.h"

typedef ptope::StackedIterator<ptope::PolytopeRebaser, ptope::PolytopeExtender,
					ptope::PolytopeCandidate> PCtoL3;
typedef ptope::CombinedCheck2<ptope::AngleCheck, true,
				ptope::DuplicateColumnCheck, false> Check;
typedef ptope::FilteredIterator<PCtoL3, ptope::PolytopeCandidate, Check, true> L3F;

typedef ptope::CombinedCheck3<ptope::AngleCheck, true,
				ptope::DuplicateColumnCheck, false, ptope::BloomPCCheck, true> Check1;
typedef ptope::CombinedCheck2<Check1, true, ptope::ParabolicCheck, false> Check2;

typedef ptope::PolytopeExtender L0toL1;
typedef ptope::FilteredIterator<L0toL1, ptope::PolytopeCandidate, Check2, true> L1F;
typedef ptope::FilteredIterator<L1F, ptope::PolytopeCandidate, ptope::PolytopeCheck, false> L1NoP;

typedef ptope::StackedIterator<L1NoP, ptope::PolytopeExtender, ptope::PolytopeCandidate> L1toL2;
typedef ptope::FilteredIterator<L1toL2, ptope::PolytopeCandidate, Check2, true> L2F;
typedef ptope::FilteredIterator<L2F, ptope::PolytopeCandidate, ptope::PolytopeCheck, false> L2NoP;

namespace ptope {
TEST(PolytopeExtender, RealCount2) {
	PolytopeCandidate p({{ 1, -.5 }, {-.5, 1} });
	Angles::get().set_angles({2,3});
	PolytopeExtender ext(p);
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(0, count);
	EXPECT_FALSE(ext.has_next());
}
TEST(PolytopeExtender, RealCount3) {
	PolytopeCandidate p({{ 1, -.5, 0 }, {-.5, 1, -.5 }, { 0, -.5, 1 } });
	Angles::get().set_angles({2,3});
	PolytopeExtender ext(p);
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(3, count);
	EXPECT_FALSE(ext.has_next());
}
TEST(PolytopeExtender, HypCount3) {
	PolytopeCandidate p({{ 1, -.5, 0 }, {-.5, 1, -.5 }, { 0, -.5, 1 } });
	PolytopeCandidate q = p.extend_by_inner_products({ -.5, -.5, -.5});
	Angles::get().set_angles({2,3});
	PolytopeExtender ext(q);
	int count = 0;
	while(ext.has_next()) {
		ext.next();
		++count;
	}
	EXPECT_EQ(3, count);
	EXPECT_FALSE(ext.has_next());
}
TEST(PolytopeExtender, NoStop) {
	PolytopeCandidate p(elliptic_factory::type_a(3));
	Angles::get().set_angles({2, 3, 4, 5, 8});
	PolytopeExtender ext(p);
	int count = 0;
	while(ext.has_next()) {
		++count;
		ext.next();
	}
	EXPECT_EQ(115, count);
	EXPECT_FALSE(ext.has_next());
	EXPECT_FALSE(ext.has_next());
	EXPECT_FALSE(ext.has_next());
	EXPECT_FALSE(ext.has_next());
}
/* Stacked Iterator assumes that the first call to next will always work, and
 * does not check has_next. This behaviour is mimicked here. */
TEST(PolytopeExtender, StackedIterMock) {
	using ptope::calc::min_cos_angle;
	Angles::get().set_angles({2, 3, 4, 5, 8, 10});
	PolytopeCandidate p(elliptic_factory::type_a(4));
	auto q = p.extend_by_inner_products({ 0, -.5, min_cos_angle(10), min_cos_angle(10) });
	auto r = q.extend_by_inner_products({ min_cos_angle(10), min_cos_angle(8), 0, 0});
	r.rebase_vectors({ 0, 1, 2, 5 });
	PolytopeExtender ext(r);
	ASSERT_TRUE(ext.next().valid());
}
TEST(ListIterator, Particular) {
	using ptope::calc::min_cos_angle;
	Angles::get().set_angles({2, 3, 4, 5, 8, 10});
	PolytopeCandidate p(elliptic_factory::type_b(4));
	auto q = p.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	auto r = q.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });

	ASSERT_TRUE(q.valid());
	ASSERT_TRUE(r.valid());

	PolytopeCandidate copy(r);
	copy.rebase_vectors({ 1, 2, 4, 5});
	PolytopeCandidate c1 = copy.extend_by_inner_products({ min_cos_angle(4), 0, 0, 0 });
	ASSERT_TRUE(c1.valid());
	PolytopeCheck pc;
	EXPECT_FALSE(pc(c1));
	DuplicateColumnCheck dc;
	EXPECT_FALSE(dc(c1));
	AngleCheck ac;
	EXPECT_TRUE(ac(c1));

	PolytopeCandidate c2 = copy.extend_by_inner_products({ 0, -.5, min_cos_angle(8), 0 });
	ASSERT_TRUE(c2.valid());

	arma::vec v1 = c1.vector_family().get(6);
	arma::vec v2 = c2.vector_family().get(6);
	double ang = calc::mink_inner_prod(v1, v2);
	EXPECT_FALSE(ang < -1);
	EXPECT_TRUE(ang < 1e-10);
	const auto & a = Angles::get().inner_products();
	comparator::DoubleLess dless;
	EXPECT_TRUE(std::binary_search(a.begin(), a.end(), ang, dless));
	EXPECT_NEAR(0, ang, 1e-10);

	L3F l3 = L3F(PCtoL3(r));
	std::set<arma::vec, comparator::VecLess> vecset;
	while(l3.has_next()) {
		auto & n = l3.next();
		auto v = n.vector_family().get(6);
		vecset.insert(v);
	}
	EXPECT_TRUE(vecset.find(v1) != vecset.end());
	EXPECT_TRUE(vecset.find(v2) != vecset.end());
}
TEST(ListIterator, Gen4Dim) {
	L0toL1 l1(elliptic_factory::type_b(4));
	L1F l1f(std::move(l1));
	L1NoP l1np(std::move(l1f));
	L1toL2 l2(std::move(l1np));
	L2F l2f(std::move(l2));
	L2NoP l2np(std::move(l2f));

	using ptope::calc::min_cos_angle;
	Angles::get().set_angles({2, 3, 4, 5, 8, 10});
	PolytopeCandidate p(elliptic_factory::type_b(4));
	auto q = p.extend_by_inner_products({ 0, min_cos_angle(8), 0, 0 });
	ASSERT_TRUE(q.valid());
	auto r = q.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	ASSERT_TRUE(r.valid());

	MEquivEqual eq;
	bool found = false;
	while(l2np.has_next()) {
		const PolytopeCandidate & n = l2np.next();
		if(eq(n.gram(), r.gram())) { 
			found = true;
			break;
		}
	}
	EXPECT_TRUE(found);
}
}

