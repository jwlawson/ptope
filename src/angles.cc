/*
 * angles.cc
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
#include "angles.h"

#include <cmath>

#include "armadillo" //Only used for pi, possibly use a different pi?

namespace ptope {
const Angles::InnerProducts &
Angles::inner_products() const {
	return _products;
}
void
Angles::set_angles(const PiSubmultiples & angles) {
	_products = angles_to_prods(angles);
}
std::vector<double>
Angles::angles_to_prods(const PiSubmultiples & angles) {
	std::vector<double> result(angles.size());
	for(std::size_t i = 0, max = angles.size(); i < max; ++i) {
		if(angles[i] == 2) {
			result[i] = 0;
		} else {
			result[i] = -std::cos(arma::datum::pi/angles[i]);
		}
	}
	std::sort(result.begin(), result.end());
	return result;
}
}

