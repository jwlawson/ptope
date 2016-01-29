/*
 * nullspace.h
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
#ifndef _PTOPE_NULLSPACE_H_
#define _PTOPE_NULLSPACE_H_

#include <armadillo>

namespace ptope {
/** 
 * Find the vector in the nullspace of A matrix with nullity 1.
 */
struct Nullspace {
public:
	bool
	operator()(arma::vec & out, const arma::mat & A);
private:
	arma::mat _qr_matrix;
	arma::podarray<double> _qr_tau;
	arma::podarray<double> _qr_work;
};
}
#endif

