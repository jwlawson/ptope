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
	auto a = angles_to_prods(angles);
	_products = std::move(a.first);
	_multiples = std::move(a.second);
}
std::pair<Angles::InnerProducts, Angles::ProdToMultiples>
Angles::angles_to_prods(const PiSubmultiples & angles) {
	std::pair<InnerProducts, ProdToMultiples> result =
		std::make_pair(InnerProducts(angles.size()), ProdToMultiples());
	for(std::size_t i = 0, max = angles.size(); i < max; ++i) {
		if(angles[i] == 2) {
			result.first[i] = 0;
			result.second[0.0] = 2;
		} else {
			result.first[i] = -std::cos(arma::datum::pi/angles[i]);
			result.second[result.first[i]] = angles[i];
		}
	}
	std::sort(result.first.begin(), result.first.end());
	return result;
}
unsigned int
Angles::inner_product(const double & d) const {
	unsigned int result;
	const auto iter = _multiples.find(d);
	if(iter != _multiples.end()) {
		result = (*iter).second;
	} else {
		result = 0;
	}
	return result;
}
}

