/*
 * compatibility_info.h
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
/**
 * Shows which vectors have compatible angles with other vectors.
 *
 * Compute the compatibilities using the from(VectorSet) method, then this set
 * will contain indices of the vectors which are compatible.
 *
 * Two vectors are considered compatible if the angle between the hyperplanes
 * normal to the vectors is one of the specified angles in Angles.get(), or the
 * hyperplanes do not intersect. (i.e. AngleCheck ( mink_inner_prod( v1, v2 ) )
 * is true ).
 */
#pragma once
#ifndef _PTOPE_COMPATIBILITY_SET_H_
#define _PTOPE_COMPATIBILITY_SET_H_

#include "angle_check.h"
#include "gram_matrix.h"
#include "vector_set.h"

#include <boost/dynamic_bitset.hpp>

namespace ptope {
class CompatibilityInfo {
public:
	CompatibilityInfo() = default;
	/**
	 * Compute which vectors in the provided VectorSet are compatible with each
	 * other.
	 */
	void from( VectorSet const& vectors );
	/**
	 * Check if the two vectors at the specified indices in the vector set are
	 * compatible.
	 */
	bool are_compatible( std::size_t const& lhs, std::size_t const& rhs) const;
	/**
	 * Get the next index after `prev` corresponding to the next vector in the
	 * vector set which is compatible to the vector corresponding to the specified
	 * `vector` index.
	 *
	 * If there is no next compatible index, then the vector index will be
	 * returned.
	 */
	std::size_t next_compatible_to( std::size_t vector, std::size_t prev ) const;

private:
	AngleCheck m_check;
	GramMatrix m_gram;
	boost::dynamic_bitset<> m_bits;
};
}
#endif
