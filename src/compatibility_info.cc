/*
 * compatibility_info.cc
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

namespace ptope {
void
CompatibilityInfo::from( VectorSet const& vectors ) {
}
bool
CompatibilityInfo::are_compatible( std::size_t const& lhs, std::size_t const& rhs) const {
	return false;
}
std::size_t
CompatibilityInfo::next_compatible_to( std::size_t vector, std::size_t prev ) const {
	return vector;
}
}

