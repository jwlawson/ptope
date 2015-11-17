/*
 * elliptic_factory_test.cc
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
#include "elliptic_factory.h"

#include <gtest/gtest.h>

namespace ptope {
namespace {
const double sqrt2 = std::sqrt(2.0);
const double min_recip_sqrt2 = -sqrt2 / 2;
const double min_cos_pi_5 = -(1+std::sqrt(5))/4;
}
TEST(EllipticFactory, Type_A) {
	arma::mat a2 = { { 1, -.5 }, { -.5, 1 } };
	arma::mat fac_a2 = elliptic_factory::type_a(2);
	arma::mat diff_2 = a2 - fac_a2;
	for(auto it = diff_2.begin(); it != diff_2.end(); ++it) {
		EXPECT_DOUBLE_EQ(0, *it);
	}

	arma::mat a3 = { { 1, -.5, 0 }, { -.5, 1, -.5 }, { 0, -.5, 1 } };
	arma::mat fac_a3 = elliptic_factory::type_a(3);
	arma::mat diff_3 = a3 - fac_a3;
	for(auto it = diff_3.begin(); it != diff_3.end(); ++it) {
		EXPECT_DOUBLE_EQ(0, *it);
	}
}
TEST(EllipticFactory, Type_B) {
	arma::mat b3 = { { 1, min_recip_sqrt2, 0 }, { min_recip_sqrt2, 1, -.5 },
		{ 0, -.5, 1 } };
	arma::mat fac_b3 = elliptic_factory::type_b(3);
	arma::mat diff_3 = b3 - fac_b3;
	for(auto it = diff_3.begin(); it != diff_3.end(); ++it) {
		EXPECT_DOUBLE_EQ(0, *it);
	}
	arma::mat b4 = { { 1, min_recip_sqrt2, 0, 0 }, { min_recip_sqrt2, 1, -.5, 0 },
		{ 0, -.5, 1, -.5 }, { 0, 0, -.5, 1 } };
	arma::mat fac_b4 = elliptic_factory::type_b(4);
	arma::mat diff_4 = b4 - fac_b4;
	for(auto it = diff_4.begin(); it != diff_4.end(); ++it) {
		EXPECT_DOUBLE_EQ(0, *it);
	}
}
TEST(EllipticFactory, Type_D) {
	arma::mat d4 = { { 1, 0, -.5, 0 }, { 0, 1, -.5, 0 }, { -.5, -.5, 1, -.5 },
		{ 0, 0, -.5, 1 } };
	arma::mat fac_d4 = elliptic_factory::type_d(4);
	auto e_it = d4.begin();
	for( auto it = fac_d4.begin(); it != fac_d4.end(); ++it) {
		EXPECT_DOUBLE_EQ(*e_it, *it);
		++e_it;
	}
}
TEST(EllipticFactory, Type_E) {
	arma::mat e6 = { { 1, -.5, 0, 0, 0, 0 }, { -.5, 1, -.5, 0, 0, 0 },
		{ 0, -.5, 1, -.5, -.5, 0 }, { 0, 0, -.5, 1, 0, 0 }, { 0, 0, -.5, 0, 1, -.5 },
		{ 0, 0, 0, 0, -.5, 1 } };
	arma::mat fac_e6 = elliptic_factory::type_e(6);
	auto e_it = e6.begin();
	for( auto it = fac_e6.begin(); it != fac_e6.end(); ++it) {
		EXPECT_DOUBLE_EQ(*e_it, *it);
		++e_it;
	}
}
TEST(EllipticFactory, Type_F) {
	arma::mat f4 = { { 1, -.5, 0, 0 }, { -.5, 1, min_recip_sqrt2, 0 },
		{ 0, min_recip_sqrt2, 1, -.5 }, { 0, 0, -.5, 1 } };
	arma::mat fac_f4 = elliptic_factory::type_f(4);
	arma::mat diff_4 = f4 - fac_f4;
	for(auto it = diff_4.begin(); it != diff_4.end(); ++it) {
		EXPECT_DOUBLE_EQ(0, *it);
	}
}
TEST(EllipticFactory, Type_G) {
	arma::mat g2_2 = { { 1, 0 }, { 0, 1 } };
	arma::mat fac_g2_2 = elliptic_factory::type_g(2,2);
	arma::mat diff_1 = g2_2 - fac_g2_2;
	for(auto it = diff_1.begin(); it != diff_1.end(); ++it) {
		EXPECT_NEAR(0, *it, 1e-15);
	}
	arma::mat g2_3 = { { 1, -.5 }, { -.5, 1 } };
	arma::mat fac_g2_3 = elliptic_factory::type_g(2,3);
	arma::mat diff_2 = g2_3 - fac_g2_3;
	for(auto it = diff_2.begin(); it != diff_2.end(); ++it) {
		EXPECT_NEAR(0, *it, 1e-15);
	}
	arma::mat g2_5 = { { 1, min_cos_pi_5 }, { min_cos_pi_5, 1 } };
	arma::mat fac_g2_5 = elliptic_factory::type_g(2,5);
	arma::mat diff_3 = g2_5 - fac_g2_5;
	for(auto it = diff_3.begin(); it != diff_3.end(); ++it) {
		EXPECT_NEAR(0, *it, 1e-15);
	}
}
TEST(EllipticFactory, Type_H) {
	arma::mat h3 = { { 1, min_cos_pi_5, 0 }, { min_cos_pi_5, 1, -.5 },
		{ 0, -.5, 1 } };
	arma::mat fac_h3 = elliptic_factory::type_h(3);
	arma::mat diff_3 = h3 - fac_h3;
	for(auto it = diff_3.begin(); it != diff_3.end(); ++it) {
		EXPECT_DOUBLE_EQ(0, *it);
	}
	arma::mat h4 = { { 1, min_cos_pi_5, 0, 0 }, { min_cos_pi_5, 1, -.5, 0 },
		{ 0, -.5, 1, -.5 }, { 0, 0, -.5, 1 } };
	arma::mat fac_h4 = elliptic_factory::type_h(4);
	arma::mat diff_4 = h4 - fac_h4;
	for(auto it = diff_4.begin(); it != diff_4.end(); ++it) {
		EXPECT_DOUBLE_EQ(0, *it);
	}
}
}

