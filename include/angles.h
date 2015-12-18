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

#include <vector>

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
	InnerProducts
	static angles_to_prods(const PiSubmultiples & angles);
	/**
	 * Remove any constructors apart form private one.
	 */
	Angles(const Angles &) = delete;
	Angles(Angles &&) = delete;
	Angles & operator=(const Angles &) = delete;
	Angles & operator=(Angles &&) = delete;
	/**
	 * Get default possible inner products.
	 */
	const InnerProducts &
	inner_products() const;
	/**
	 * Set the possible default angles.
	 */
	void
	set_angles(const PiSubmultiples & angles);
private:
	Angles() : _products(angles_to_prods({ 2, 3, 4, 5, 8})) {}
	InnerProducts _products;
};
}
#endif

