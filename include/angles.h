/*
 * angles.h
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
#pragma once
#ifndef PTOPE_ANGLES_H_
#define PTOPE_ANGLES_H_

#include <map>
#include <vector>

#include "comparator.h"

namespace ptope {
/**
 * Don't know if this is the best way of providing config.
 *
 * Provide central lookup for values of inner products to use when extending
 * or checking polytope gram matrices.
 */
class Angles {
public:
	typedef std::vector<double> InnerProducts;
	typedef std::vector<unsigned int> PiSubmultiples;
	typedef std::map<double, unsigned int, comparator::DoubleLess> ProdToMultiples;
	/**
	 * Get singleton instance.
	 */
	static 
	Angles &
	get() {
		static Angles instance;
		return instance;
	}
	/**
	 * Compute inner products from angle submultiples.
	 */
	std::pair<InnerProducts, ProdToMultiples>
	static angles_to_prods(const PiSubmultiples & angles);
	/**
	 * Remove any constructors apart form private one.
	 */
	Angles(const Angles &) = delete;
	Angles(Angles &&) = delete;
	Angles & operator=(const Angles &) = delete;
	Angles & operator=(Angles &&) = delete;
	/**
	 * Get default possible inner products. This has been sorted, so doesn't need
	 * to be sorted again.
	 */
	const InnerProducts &
	inner_products() const;
	/**
	 * Set the possible default angles.
	 */
	void
	set_angles(const PiSubmultiples & angles);
	/**
	 * Check if given inner product is included in list of products and if so
	 * return the corresponding angle submultiple.
	 */
	unsigned int
	inner_product(const double & d) const;
private:
	Angles() 
	: _products(angles_to_prods({ 2, 3, 4, 5, 8}).first),
		_multiples(angles_to_prods({2, 3, 4, 5, 8}).second) {}
	InnerProducts _products;
	ProdToMultiples _multiples;
};
}
#endif

