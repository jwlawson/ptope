/*
 * gram_matrix.h
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
/**
 */
#pragma once
#ifndef PTOPE_GRAM_MATRIX_H_
#define PTOPE_GRAM_MATRIX_H_

#include "vector_set.h"

#include <tuple>

#include <boost/iterator/iterator_facade.hpp>

/**
 * A Gram matrix for a set of vectors.
 *
 * Under the hood the matrix is stored in RFP form - that is rectangular full
 * packed ( see http://www.netlib.org/lapack/lawnspdf/lawn199.pdf ). The matrix
 * is stored in normal (not transposed) lower format.
 *
 * The class provides some human accessible methods for getting inner products
 * from the matrix, as well as some more direct access to the RFP matrix for
 * anyone who knows what they are doing.
 */
namespace ptope {
class GramMatrix {
	static const char L;
	static const char N;
	static const char T;
	static const double one;
	static const double minus_one;
	static const arma::blas_int int_one;
public:
	class const_iterator;
	struct Entry {
		std::size_t row;
		std::size_t col;
		double val;
	};
	GramMatrix() = default;
	/**
	 * Compute the gram matrix of the provided vectors.
	 */
	void from( VectorSet const& vectors );
	/**
	 * Get the row, column and inner product at the specified RFP index.
	 */
	Entry at_rfp( std::size_t const& index ) const;
	/**
	 * Get the inner product at the (row, col) position in the matrix.
	 *
	 * As the matrix is symmetric, at( row , col ) == at( col , row ).
	 */
	double at( std::size_t const& row, std::size_t const& col ) const;
	/** Get just the entry at the specified index. */
	double raw_rfp( std::size_t const& index ) const;

	const_iterator begin() const;
	const_iterator end() const;
private:
	std::size_t m_nvecs;
	std::size_t m_nelems;
	std::size_t m_rfp_ncols;
	std::size_t m_rfp_nrows;
	arma::podarray<double> m_matrix;

	/** Get the smallest size needed to hold all inner products of vectors. */
	std::size_t
	priv_min_product_size( VectorSet const& vectors );
	/** Get the row, col pair for an rfp index. */
	std::pair<std::size_t, std::size_t>
	priv_rfp_index_to_ij( std::size_t const& index ) const;
	std::size_t
	priv_ij_to_rfp_index( std::size_t const& row, std::size_t const& col ) const;
	/** Prepare to compute products.
	 * Make sure the workspace is large enough etc.
	 */
	void
	priv_products_prepare( ptope::VectorSet const& vectors );
	/** Compute the inner products of the vectors. */
	void
	priv_products_compute( ptope::VectorSet const& vectors );
};

inline
void
GramMatrix::from( ptope::VectorSet const& vectors ) {
	priv_products_prepare( vectors );
	priv_products_compute( vectors );
}
inline
GramMatrix::Entry
GramMatrix::at_rfp( std::size_t const& index ) const {
	Entry result;
	std::tie(result.row, result.col) = priv_rfp_index_to_ij( index );
	result.val = m_matrix[index];
	return result;
}
inline
double
GramMatrix::at( std::size_t const& row, std::size_t const& col ) const {
	return m_matrix[ priv_ij_to_rfp_index( row, col ) ];
}
inline
std::size_t
GramMatrix::priv_min_product_size( ptope::VectorSet const& vectors ) {
	std::size_t num = vectors.size();
	return ( num * ( num + 1 ) ) / 2;
}

class GramMatrix::const_iterator 
	: public boost::iterator_facade< 
			  GramMatrix::const_iterator
			, GramMatrix::Entry const
			, std::random_access_iterator_tag > {
	typedef const_iterator self_type;
public:
	const_iterator() = default;
	const_iterator( self_type const& ) = default;
	const_iterator( self_type && ) = default;
	const_iterator( std::size_t index, GramMatrix const& gram )
		: m_index{ index }
		, m_entry{}
		, m_gram(gram)
	{}

private:
	friend class boost::iterator_core_access;

	std::size_t m_index;
	GramMatrix::Entry m_entry;
	GramMatrix const& m_gram;

	reference dereference()
	{ m_entry = m_gram.at_rfp(m_index); return m_entry; }
	bool equal( self_type const& other )
	{ return m_index == other.m_index; }
	void increment()
	{ ++m_index; }
	void decrement()
	{ --m_index; }
	void advance( difference_type offset )
	{ m_index += offset; }
	difference_type distance_to( self_type const& other )
	{ return other.m_index - m_index; }
};
#ifdef ARMA_USE_BLAS
#ifndef ARMA_BLAS_CAPITALS
	#define arma_dsfrk dsfrk
#else
	#define arma_dsfrk DSFRK
#endif
namespace blas {
extern "C" {
	/**
	 * Level 3 BLAS like routine for C in RFP Format.
	 *
	 * DSFRK performs one of the symmetric rank--k operations
	 *
	 *		C := alpha*A*A**T + beta*C,
	 *
	 * or
	 *
	 *		C := alpha*A**T*A + beta*C,
	 *
	 * where alpha and beta are real scalars, C is an n--by--n symmetric
	 * matrix and A is an n--by--k matrix in the first case and a k--by--n
	 * matrix in the second case.
	 */
	void arma_fortran(arma_dsfrk)(
		const char * TRANSR,
		const char * UPLO,
		const char * TRANS,
		const arma::blas_int * N,
		const arma::blas_int * K,
		const double * ALPHA,
		const double * A,
		const arma::blas_int * LDA,
		const double * BETA,
		double * C
	);
}
}
#endif
}
#endif

