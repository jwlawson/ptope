/*
 * polytope_candidate_test.cc
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
#include "polytope_candidate.h"
#include "elliptic_factory.h"

#include <gtest/gtest.h>

namespace ptope {
namespace {
constexpr double error = 1e-15;
double min_cos_angle(uint mult) {
	return -std::cos(arma::datum::pi/mult);
}
}
TEST(PolytopeCandidate, FirstExtend) {
	PolytopeCandidate p({ { 1, -.5 }, { -.5, 1 } });
	auto q = p.extend_by_inner_products({ min_cos_angle(4), min_cos_angle(3) });
	ASSERT_TRUE(q.valid());
	arma::mat exp = 
	{ { 1, -.5, min_cos_angle(4) },
		{ -.5, 1, min_cos_angle(3) },
		{ min_cos_angle(4), min_cos_angle(3), 1 } };
	arma::mat diff = q.gram() - exp;
	for(const double & val : diff) {
		EXPECT_NEAR(0.0, val, error);
	}
}
TEST(PolytopeCandidate, EuclideanInvalid) {
	PolytopeCandidate p({ { 1, -.5 }, { -.5, 1 } });
	auto q = p.extend_by_inner_products({ -.5, -.5 });
	EXPECT_FALSE(q.valid());
}
TEST(PolytopeCandidate, EsselmannExample) {
	PolytopeCandidate p({ { 1, -.5, 0, 0 }, 
												{ -.5, 1, min_cos_angle(4), 0 }, 
												{ 0, min_cos_angle(4), 1, -.5 }, 
												{ 0, 0, -.5, 1 } });
	auto q = p.extend_by_inner_products({ 0, 0, 0, min_cos_angle(8) });
	ASSERT_TRUE(q.valid());
	auto r = q.extend_by_inner_products({ min_cos_angle(8), 0, 0, 0 });
	ASSERT_TRUE(r.valid());
	arma::mat exp = 
	{ {  1, -.5, 0, 0, 0, min_cos_angle(8) }, 
		{ -.5, 1, min_cos_angle(4), 0, 0, 0 }, 
		{ 0, min_cos_angle(4), 1, -.5, 0, 0 }, 
		{ 0, 0, -.5, 1, min_cos_angle(8), 0 },
		{ 0, 0, 0, min_cos_angle(8), 1, 0 },
		{ min_cos_angle(8), 0, 0, 0, 0, 1 } };
	arma::mat diff = r.gram() - exp;
	for(const double & val : diff) {
		EXPECT_NEAR(0.0, val, error);
	}
}
TEST(PolytopeCandidate, LannerExample) {
	PolytopeCandidate p({ { 1, -.5, 0 }, { -.5, 1, min_cos_angle(5) },
			{ 0, min_cos_angle(5), 1 }});
		auto q = p.extend_by_inner_products({ 0, 0, -.5});
		ASSERT_TRUE(q.valid());
		arma::mat exp = { { 1, -.5, 0, 0 }, { -.5, 1, min_cos_angle(5), 0},
			{ 0, min_cos_angle(5), 1, -.5 }, { 0, 0, -.5, 1} };
	arma::mat diff = q.gram() - exp;
	for(const double & val : diff) {
		EXPECT_NEAR(0.0, val, error);
	}
}
TEST(PolytopeCandidate, StrangeNaN) {
	arma::mat a6 = elliptic_factory::type_a(6);
	PolytopeCandidate p(a6);
	PolytopeCandidate q = p.extend_by_inner_products({ min_cos_angle(5), 0, 0, 0, 0, 0 } );
	PolytopeCandidate r = q.extend_by_inner_products({ 0, 0, 0, 0, min_cos_angle(4), min_cos_angle(3) });
	EXPECT_FALSE(q.gram().has_nan());
}
TEST(PolytopeCandidate, SaveLoad) {
	PolytopeCandidate p({ { 1, 0 }, { 0, 1 } });
	p.extend_by_inner_products({ -.5, -.5 });
	std::stringstream ss;
	p.save(ss);
	
	PolytopeCandidate q;
	q.load(ss);

	const arma::mat & p_g = p.gram();
	const arma::mat & q_g = q.gram();
	for(arma::uword i = 0; i < p_g.size(); ++i) {
		EXPECT_DOUBLE_EQ(p_g(i), q_g(i));
	}
	EXPECT_EQ(p.valid(), q.valid());
}
}

