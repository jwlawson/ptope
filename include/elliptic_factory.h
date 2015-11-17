/*
 * elliptic_factory.h
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
 * Set of factory methods to construct different elliptic coxeter matrices.
 */
#pragma once
#ifndef PTOPE_ELLIPTIC_FACTOR_H_
#define PTOPE_ELLIPTIC_FACTOR_H_

#include <armadillo>

namespace ptope {
namespace elliptic_factory {
/**
 * Return the Gram matrix of Coxeter diagram of type An, with specified number
 * of vertices.
 */
arma::mat type_a(const uint size);
/**
 * Return the Gram matrix of Coxeter diagram of type Bn = Cn, with specified number
 * of vertices.
 */
arma::mat type_b(const uint size);
/**
 * Return the Gram matrix of Coxeter diagram of type Dn, with specified number
 * of vertices.
 */
arma::mat type_d(const uint size);
/**
 * Return the Gram matrix of Coxeter diagram of type En, with specified number
 * of vertices.
 *
 * Only valid for sizes 6, 7 and 8. Otherwise a 0x0 matrix is returned.
 */
arma::mat type_e(const uint size);
/**
 * Return the Gram matrix of Coxeter diagram of type Fn, with specified number
 * of vertices.
 *
 * Only valid for size 4. Otherwise a 0x0 matrix is returned.
 */
arma::mat type_f(const uint size);
/**
 * Return the Gram matrix of Coxeter diagram of type Gn, with specified number
 * of vertices, and the label on the edge.
 *
 * Only valid for size 2. Otherwise a 0x0 matrix is returned.
 */
arma::mat type_g(const uint size, const uint label);
/**
 * Return the Gram matrix of Coxeter diagram of type Hn, with specified number
 * of vertices.
 *
 * Only valid for sizes 3 and 4. Otherwise a 0x0 matrix is returned.
 */
arma::mat type_h(const uint size);
}
}
#endif

