/*
 * elliptic_generator.h
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
/**
 * Class to generate all elliptic Coxeter schemes for a given size.
 *
 * Uses a java iterator style interface to generate each scheme sequentially.
 */
#pragma once
#ifndef PTOPE_ELLIPTIC_GENERATOR_H_
#define PTOPE_ELLIPTIC_GENERATOR_H_

#include "elliptic_factory.h"

namespace ptope {
class EllipticGenerator {
	public:
		/**
		 * Create an instance with the specified size.
		 */
		EllipticGenerator(uint size);
		/**
		 * Check whether there are further matrices to generate. Returns true if a
		 * subsequent call to next() will return a valid matrix.
		 */
		bool hasNext();
		/**
		 * Generates and returns the next elliptic matrix.
		 */
		arma::mat next();
}
#endif

