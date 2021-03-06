/*
 * gram_matrix.cc
 * Copyright 2016 John Lawson
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
#include "gram_matrix.h"


namespace ptope {
const char GramMatrix::L = 'L';
const char GramMatrix::N = 'N';
const char GramMatrix::T = 'T';
const double GramMatrix::one = 1.0;
const double GramMatrix::minus_one = -1.0;
const arma::blas_int GramMatrix::int_one = 1;

std::pair<std::size_t, std::size_t>
GramMatrix::priv_rfp_index_to_ij( std::size_t const& index ) const {
	std::size_t const rfp_col = index / m_rfp_nrows;
	bool const is_even = m_nvecs % 2 == 0;
	std::size_t row_cutoff = rfp_col;
	if( is_even ) { ++row_cutoff; }
	std::size_t const rfp_row = index % m_rfp_nrows;
	std::size_t row;
	std::size_t col;
	if( rfp_row < row_cutoff ) {
		row = m_rfp_ncols + rfp_col;
		col = m_rfp_ncols + rfp_row;
		if( !is_even ) { --row; }
	} else {
		col = rfp_col;
		row = rfp_row;
		if( is_even ) { --row; }
	}
	return { row, col };
}
std::size_t
GramMatrix::priv_ij_to_rfp_index( std::size_t const& row ,
		std::size_t const& col ) const {
	std::size_t result;
	if( row < col ) {
		result = priv_ij_to_rfp_index( col, row );
	} else {
		std::size_t rfp_row;
		std::size_t rfp_col;
		std::size_t const col_threshold = m_rfp_ncols;
		if( col < col_threshold ) {
			rfp_row = row;
			if( m_nvecs % 2 == 0 ) { ++rfp_row; }
			rfp_col = col;
		} else {
			rfp_row = col - col_threshold;
			rfp_col = row - col_threshold;
			if( m_nvecs % 2  == 1 ) { ++rfp_col; }
		}
		result = rfp_col * m_rfp_nrows + rfp_row;
	}
	return result;
}
void
GramMatrix::priv_products_prepare( VectorSet<double> const& vectors ) {
	m_nvecs = vectors.size();
	m_nelems = priv_min_product_size( vectors );
	m_rfp_ncols = ( m_nvecs  + 1 ) / 2;
	m_rfp_nrows = m_nvecs  + 1 - ( m_nvecs  % 2 );
	m_matrix.set_min_size( m_nelems );
	m_matrix.zeros();
}
void
GramMatrix::priv_products_compute( VectorSet<double> const& vectors ) {
	std::size_t const total_dim = vectors.dimension();
	std::size_t const real_dim = total_dim - 1;
	double const * vector_ptr = vectors.memptr();
	double const * hyp_part_ptr = vector_ptr + real_dim;
	arma::blas_int n = m_nvecs;
	arma::blas_int k = real_dim;
	arma::blas_int lda = total_dim;
	double * out_ptr = m_matrix.memptr();

	/* Compute hyperbolic part of product */
	blas::arma_fortran(arma_dsfrk)(&N, &L, &T, &n, &int_one, &one, hyp_part_ptr, &lda,
			&one, out_ptr);
	/* Compute inner product of euclidean part and subtract the hyperbolic
	 * product. */
	blas::arma_fortran(arma_dsfrk)(&N, &L, &T, &n, &k, &one, vector_ptr, &lda,
			&minus_one, out_ptr);
}
}

