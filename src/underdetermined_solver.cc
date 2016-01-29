/*
 * underdetermined_solver.cc
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
#include "underdetermined_solver.h"

namespace ptope {
/* Simplification of arma::solve */
bool
UDSolver::operator()(arma::vec & out, const arma::mat & A, const arma::vec & B) {
	using arma::uword;
	using arma::blas_int;
	_a = A;
	const uword A_n_rows = A.n_rows;
	const uword A_n_cols = A.n_cols;

	const uword B_n_rows = B.n_rows;
	const uword B_n_cols = B.n_cols;

	char trans = 'N';
	blas_int  m     = blas_int(A_n_rows);
	blas_int  n     = blas_int(A_n_cols);
	blas_int  lda   = blas_int(A_n_rows);
	blas_int  ldb   = blas_int(A_n_cols);
	blas_int  nrhs  = blas_int(B_n_cols);
	blas_int  lwork = 3 * ( (std::max)(blas_int(1), m + (std::max)(m,nrhs)) );
	blas_int  info  = 0;

	out.set_size(A_n_cols, B_n_cols);

	double* sol_ptr = out.memptr();
	arma::arrayops::copy( sol_ptr, B.memptr(), B_n_rows );
	for(uword row=B_n_rows; row<A_n_cols; ++row) {
		sol_ptr[row] = double(0);
	}
	_work.set_min_size(static_cast<uword>(lwork));

	arma::lapack::gels<double>( &trans, &m, &n, &nrhs, _a.memptr(), &lda,
			out.memptr(), &ldb, _work.memptr(), &lwork, &info );

	return (info == 0);
}
}

