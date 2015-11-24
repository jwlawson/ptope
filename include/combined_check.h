/*
 * combined_check.h
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
#ifndef PTOPE_COMBINED_CHECK_H_
#define PTOPE_COMBINED_CHECK_H_

#include "polytope_candidate.h"

namespace ptope {
template <class Chk1, bool pos1, class Chk2, bool pos2>
class CombinedCheck2 {
public:
	bool operator()(const PolytopeCandidate & p) {
		return (_chk1(p) == pos1) && (_chk2(p) == pos2);
	}
private:
	Chk1 _chk1;
	Chk2 _chk2;
};
template <class Chk1, bool pos1, class Chk2, bool pos2, class Chk3, bool pos3>
class CombinedCheck3 {
public:
	bool operator()(const PolytopeCandidate & p) {
		return (_chk1(p) == pos1) && (_chk2(p) == pos2) && (_chk3(p) == pos3);
	}
private:
	Chk1 _chk1;
	Chk2 _chk2;
	Chk3 _chk3;
};
}
#endif

