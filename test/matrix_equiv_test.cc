#include "matrix_equiv.h"

#include <gtest/gtest.h>

#include "elliptic_factory.h"
#include "polytope_extender.h"

namespace ptope {
namespace {
double min_cos_angle(uint mult) {
	return -std::cos(arma::datum::pi/mult);
}
}
TEST(MatrixEquiv, Equal) {
	MEquivEqual eq;
	arma::mat a = { { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 } };
	arma::mat b = { { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 } };
	EXPECT_TRUE(eq(a,b));
	b.swap_cols(0,1);
	b.swap_rows(0,1);
	EXPECT_TRUE(eq(a,b));
	b.swap_cols(0,2);
	b.swap_rows(0,2);
	EXPECT_TRUE(eq(a,b));
}
TEST(MatrixEquiv, NEqual) {
	MEquivEqual eq;
	arma::mat a = { { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 } };
	arma::mat b = { { 9, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 } };
	EXPECT_FALSE(eq(a,b));
	b.swap_cols(0,1);
	b.swap_rows(0,1);
	EXPECT_FALSE(eq(a,b));
	b.swap_cols(0,2);
	b.swap_rows(0,2);
	EXPECT_FALSE(eq(a,b));
}
TEST(MatrixHash, Equal) {
	MEquivHash hash;
	arma::mat a = { { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 } };
	arma::mat b = { { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 } };
	EXPECT_EQ(hash(a), hash(b));
	b.swap_cols(0,1);
	b.swap_rows(0,1);
	EXPECT_EQ(hash(a), hash(b));
	b.swap_cols(0,2);
	b.swap_rows(0,2);
	EXPECT_EQ(hash(a), hash(b));
}
TEST(MatrixHash, NEqual) {
	MEquivHash hash;
	arma::mat a = { { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 } };
	arma::mat b = { { 9, 2, 3 }, { 2, 3, 4 }, { 3, 4, 5 } };
	EXPECT_NE(hash(a), hash(b));
	b.swap_cols(0,1);
	b.swap_rows(0,1);
	EXPECT_NE(hash(a), hash(b));
	b.swap_cols(0,2);
	b.swap_rows(0,2);
	EXPECT_NE(hash(a), hash(b));
}
TEST(MatrixEqual, Polytopes) {
	MEquivEqual e;
	MEquivHash h;
	PolytopeCandidate p(elliptic_factory::type_a(4));
	auto q = p.extend_by_inner_products({ 0, 0, 0, min_cos_angle(5) });
	ASSERT_TRUE(q.valid());

	auto s = p.extend_by_inner_products({ 0, 0, 0, min_cos_angle(5) });
	ASSERT_TRUE(s.valid());
	EXPECT_EQ(h(q.gram()), h(s.gram()));
	EXPECT_TRUE(e(q.gram(), s.gram()));

	auto r = p.extend_by_inner_products({ min_cos_angle(5), 0, 0, 0 });
	ASSERT_TRUE(r.valid());
	EXPECT_EQ(h(q.gram()), h(r.gram()));
	EXPECT_TRUE(e(q.gram(), r.gram()));
}
TEST(MatrixEqual, LargeNumber) {
	MEquivEqual e;
	MEquivHash h;
	PolytopeCandidate p(elliptic_factory::type_a(4));
	auto q = p.extend_by_inner_products({ 0, 0, 0, min_cos_angle(5) });
	ASSERT_TRUE(q.valid());
	PolytopeExtender ext(p);
	int e_count = 0;
	int h_count = 0;
	while(ext.has_next()) {
		const auto & n = ext.next();
		if(h(q.gram()) == h(n.gram())) ++h_count;
		if(e(q.gram(), n.gram())) ++e_count;
	}
	EXPECT_EQ(2, h_count);
	EXPECT_EQ(2, e_count);
}
}

