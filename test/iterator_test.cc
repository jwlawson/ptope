#include "elliptic_factory.h"
#include "filtered_iterator.h"
#include "number_dotted_check.h"
#include "polytope_candidate.h"
#include "polytope_extender.h"

#include <gtest/gtest.h>

namespace ptope {
TEST(FIterator, One) {
	PolytopeCandidate p(elliptic_factory::type_a(3));
	PolytopeExtender ext(p);
	FilteredIterator<PolytopeExtender, PolytopeCandidate, NumberDottedCheck<0>, false> f(std::move(ext));
	int count = 0;
	while(f.has_next() && count < 100) {
		++count;
		f.next();
	}
	EXPECT_EQ(0, count);
}
}

