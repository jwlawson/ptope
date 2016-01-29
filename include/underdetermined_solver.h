/*
 * underdetermined_solver.h
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
#ifndef _PTOPE_UNDERDETERMINED_SOLVER_H_
#define _PTOPE_UNDERDETERMINED_SOLVER_H_

#include <armadillo>

namespace ptope {
/** 
 * Solve an underdetermined system of linear equations, using BLAS/LAPACK dgels.
 */
struct UDSolver {
public:
	bool
	operator()(arma::vec & out, const arma::mat & A, const arma::vec & B);
private:
	arma::mat _a;
	arma::podarray<double> _work;
};
}
#endif

