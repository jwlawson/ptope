/*
 * nullspace.cc
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
#include "nullspace.h"

namespace ptope {
bool
Nullspace::operator()(arma::vec & out, const arma::mat & X) {
	using arma::uword;
	using arma::blas_int;
	/* This makes QR too big for just R, but the right size for Q. qeqrf fills QR
	 * with R in upper triangle and Q encoded in the bottom. orgqr decodes Q into
	 * the fill Q matrix. By Using the same matrix for everything we avoid having
	 * to copy all the values between Q and R between the function calls. */
	_qr_matrix.set_size(X.n_cols, X.n_cols);
	_qr_matrix.head_cols(X.n_rows) = X.t();

	const uword R_n_rows = _qr_matrix.n_rows;
	const uword R_n_cols = _qr_matrix.n_cols - 1;

	blas_int m         = static_cast<blas_int>(R_n_rows);
	blas_int n         = static_cast<blas_int>(R_n_cols);
	blas_int lwork     = 0;
	// take into account requirements of geqrf() _and_ orgqr()/ungqr()
	blas_int lwork_min = (std::max)(blas_int(1), (std::max)(m,n));
	blas_int k         = (std::min)(m,n);
	blas_int info      = 0;

	_qr_tau.set_min_size( static_cast<uword>(k) );

	double work_query[2];
	blas_int lwork_query = -1;

	arma::lapack::geqrf(&m, &n, _qr_matrix.memptr(), &m, _qr_tau.memptr(),
			&work_query[0], &lwork_query, &info);

	if(info != 0)  { return false; }

	blas_int lwork_proposed = static_cast<blas_int>(
			arma::access::tmp_real(work_query[0]) );
	lwork = (std::max)(lwork_proposed, lwork_min);
	_qr_work.set_min_size( static_cast<uword>(lwork) );

	arma::lapack::geqrf(&m, &n, _qr_matrix.memptr(), &m, _qr_tau.memptr(),
			_qr_work.memptr(), &lwork, &info);

	if(info != 0)  { return false; }

	arma::lapack::orgqr(&m, &m, &k, _qr_matrix.memptr(), &m, _qr_tau.memptr(),
			_qr_work.memptr(), &lwork, &info);

	out.set_size(R_n_rows);
	arma::arrayops::copy(out.memptr(), _qr_matrix.colptr(R_n_rows - 1), R_n_rows);
	return (info == 0);
}
}

