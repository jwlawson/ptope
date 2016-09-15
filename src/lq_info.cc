/*
 * lq_info.cc
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
#include "lq_info.h"

namespace ptope {
namespace detail {
namespace {
arma::mat __li_cache;
arma::mat __orthog_cache;
arma::podarray<double> __work;
arma::podarray<double> __tau;
arma::vec __nullvec;
}
std::unique_ptr<LQInfo>
LQInfo::compute(arma::mat const & A) {
	using arma::uword;
	using arma::blas_int;
	const uword A_n_rows = A.n_rows;
	const uword A_n_cols = A.n_cols;
	std::unique_ptr<LQInfo> result(new LQInfo(A_n_rows));
	__orthog_cache.set_size(A_n_cols, A_n_cols);
	__orthog_cache.submat(0, 0, A_n_rows - 1, A_n_cols - 1) = A;

	blas_int m(A_n_rows);
	blas_int n(A_n_cols);
	/* _lq_cache has an extra row to A, so lda = m + 1 = n */
	blas_int lda(A_n_cols);
	blas_int info(0);
	blas_int k = std::min(m,n);
	blas_int lwork_min = std::max(blas_int(1), std::max(m,n));
	double work_query[2];
	blas_int lwork_query(-1);

	__tau.set_min_size(static_cast<uword>(k));

	/* Find the optimum work space size for computing LQ */
	arma_fortran(dgelqf)(&m, &n, __orthog_cache.memptr(), &lda, __tau.memptr(),
			&work_query[0], &lwork_query, &info);

	blas_int lwork_proposed = static_cast<blas_int>(
			arma::access::tmp_real(work_query[0]) );
	blas_int lwork = (std::max)(lwork_proposed, lwork_min);
	__work.set_min_size( static_cast<uword>(lwork) );

	/* Compute LQ decomposition of A */
	arma_fortran(dgelqf)(&m, &n, __orthog_cache.memptr(), &lda, __tau.memptr(),
			__work.memptr(), &lwork, &info);
	__li_cache = arma::inv(arma::trimatl(__orthog_cache.submat(0, 0, A_n_rows -
					1, A_n_rows - 1)));
	/* Compute orthogonal matrix from LQ decomp to get nullspace */
	arma_fortran(dorglq)(&n, &n, &k, __orthog_cache.memptr(), &n, __tau.memptr(),
			__work.memptr(), &lwork, &info);

	/* Write nullspace to output. */
	result->_nullspace = __orthog_cache.row(A_n_rows).t();
	result->_qtli = __orthog_cache.head_rows(A_n_rows).t() * __li_cache;
	
	return result;
}
LQInfo::LQInfo(arma::uword const dim)
	: _qtli(dim + 1, dim),
		_nullspace(dim + 1) {}
}
}

