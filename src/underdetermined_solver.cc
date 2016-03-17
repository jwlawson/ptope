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

#include "nullspace.h"

namespace ptope {
/* Simplification of arma::solve */
bool
UDSolver::operator()(arma::vec & out, arma::vec & nullvec,
		detail::LQInfo & lqi, const arma::vec & B) {
	using arma::uword;
	using arma::blas_int;
	const arma::mat & A = lqi.lq();
	const uword A_n_rows = A.n_rows;
	const uword A_n_cols = A.n_cols;
	const uword B_n_rows = B.n_rows;

	char ntrans = 'N';
	char low = 'L';
	char trans = 'T';
	blas_int m(A_n_rows);
	blas_int n(A_n_cols);
	blas_int lda(A_n_rows);
	blas_int ldb(A_n_cols);
	blas_int nrhs(1);
	blas_int info(0);
	blas_int lwork_min = std::max(blas_int(1), std::max(m,n));
	double work_query[2];
	blas_int lwork_query(-1);

	/* Initialize out */
	out.set_size(A_n_cols);
	double* sol_ptr = out.memptr();
	arma::arrayops::copy( sol_ptr, B.memptr(), B_n_rows );
	for(uword row=B_n_rows; row<A_n_cols; ++row) {
		sol_ptr[row] = double(0);
	}

	/* Find the optimum work space size for computing solution */
	arma_fortran(arma_dormlq)(&low, &trans, &n, &nrhs, &m, lqi.lq_ptr(), &lda,
			lqi.tau_ptr(), out.memptr(), &ldb, &work_query[0], &lwork_query, &info);
	if(info != 0) return false;

	blas_int lwork_proposed = static_cast<blas_int>(
			arma::access::tmp_real(work_query[0]) );
	blas_int lwork = (std::max)(lwork_proposed, lwork_min);
	_work.set_min_size( static_cast<uword>(lwork) );

	/* Solve triangular system of equations */
	arma::dtrtrs_(&low, &ntrans, &ntrans, &m, &nrhs, lqi.lq_ptr(), &lda,
				out.memptr(), &ldb, &info);
	if(info != 0) return false;
	for(int i = m; i < n; ++i) {
		out(i) = 0;
	}
	/* Change basis usig orthogonal part of A */
	arma_fortran(arma_dormlq)(&low, &trans, &n, &nrhs, &m, lqi.lq_ptr(), &lda,
			lqi.tau_ptr(), out.memptr(), &ldb, _work.memptr(), &lwork, &info);
	if(info != 0) return false;

	/* Write nullspace to output. */
	nullvec = lqi.null();

	return (info == 0);
}
}

