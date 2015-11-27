/*
 * angle_check.h
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
#ifndef PTOPE_ANGLE_CHECK_H_
#define PTOPE_ANGLE_CHECK_H_

#include "polytope_candidate.h"

namespace ptope {
/**
 * Checks whether all entries in the last column of the matrix/polytope are one
 * of the provided angles. Returns true if all angles are valid, otherwise
 * returns false.
 */
class AngleCheck {
static const std::vector<uint> __default_mults;
static const std::vector<double> __default_inner;
public:
	AngleCheck();
	AngleCheck(const std::vector<uint> & angles);
	AngleCheck(std::vector<uint> && angles);
	bool operator()(const PolytopeCandidate & p);
	bool operator()(const arma::mat & m);
private:
	std::vector<double> _values;
};
}
#endif

