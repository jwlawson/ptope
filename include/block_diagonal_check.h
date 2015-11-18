/*
 * block_diagonal_check.h
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
 * Functor which checks whether a given matrix is block diagonal or not (after
 * some permutation). It returns true if the matrix is block diagonal.
 */
#pragma once
#ifndef PTOPE_PARABOLIC_FILTER_H_
#define PTOPE_PARABOLIC_FILTER_H_

#include <armadillo>
#include <queue>

namespace ptope {
class BlockDiagonalCheck {
	public:
		BlockDiagonalCheck() {}
		bool operator()(const arma::mat m);
	private:
		std::vector<arma::uword> _unvisited;
		std::vector<std::vector<arma::uword>> _components;
		std::queue<arma::uword> _queue;
};
}
#endif

