/*
 * compatibility_info.cc
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
#include "compatibility_info.h"

namespace ptope {
void
CompatibilityInfo::from( VectorSet<double> const& vectors ) {
	m_num_vecs = vectors.size();
	m_gram.from( vectors );
	std::size_t const num_entries = ( m_num_vecs * (m_num_vecs + 1) ) / 2;
	m_bits.resize( num_entries );
	if( m_num_vecs % 2 == 0 ) {
		priv_compatibility_even();
	} else {
		priv_compatibility_odd();
	}
}
bool
CompatibilityInfo::are_compatible( std::size_t const& lhs, std::size_t const& rhs) const {
	return m_bits.test( priv_index_from_ij( lhs , rhs ) );
}
std::size_t
CompatibilityInfo::next_compatible_to( std::size_t vector, std::size_t prev ) const {
	if( prev < vector ) { prev = vector; }
	std::size_t const base = priv_index_from_ij( vector, vector );
	std::size_t const based_prev = base + ( prev - vector );
	std::size_t const next_base = m_bits.find_next( based_prev );
	std::size_t const max_index = base + m_num_vecs - vector - 1;
	std::size_t next_ind;
	if( next_base > max_index ) {
		next_ind = vector;
	} else {
		next_ind = vector + ( next_base - base );
	}
	return next_ind;
}
void
CompatibilityInfo::priv_compatibility_even() {
	std::size_t const gram_ncols = m_num_vecs / 2;
	std::size_t const gram_nrows = m_num_vecs + 1;
	std::size_t bits_ind = 0;
	std::size_t gram_ind = 0;
	for( std::size_t col = 0; col < gram_ncols; ++col ) {
		gram_ind += col + 1;
		for( std::size_t row = 1 + col; row < gram_nrows; ++row, ++bits_ind, ++gram_ind ) {
			double const val = m_gram.raw_rfp( gram_ind );
			m_bits.set( bits_ind , m_check( val ) );
		}
	}
	for( std::size_t row = 0; row < gram_ncols; ++row ) {
		for( std::size_t col = row; col < gram_ncols; ++col, ++bits_ind ) {
			gram_ind = row + gram_nrows * col;
			double const val = m_gram.raw_rfp( gram_ind );
			m_bits.set( bits_ind , m_check( val ) );
		}
	}
}
void
CompatibilityInfo::priv_compatibility_odd() {
	std::size_t const gram_ncols = m_num_vecs / 2 + 1;
	std::size_t const gram_nrows = m_num_vecs;
	std::size_t bits_ind = 0;
	std::size_t gram_ind = 0;
	for( std::size_t col = 0; col < gram_ncols; ++col ) {
		gram_ind += col;
		for( std::size_t row = col; row < gram_nrows; ++row, ++bits_ind, ++gram_ind ) {
			double const val = m_gram.raw_rfp( gram_ind );
			m_bits.set( bits_ind , m_check( val ) );
		}
	}
	std::size_t const sub_rfp_nrows = gram_nrows - gram_ncols;
	for( std::size_t row = 0; row < sub_rfp_nrows; ++row ) {
		for( std::size_t col = 1 + row; col < gram_ncols; ++col, ++bits_ind ) {
			gram_ind = row + gram_nrows * col;
			double const val = m_gram.raw_rfp( gram_ind );
			m_bits.set( bits_ind , m_check( val ) );
		}
	}
}
std::size_t
CompatibilityInfo::priv_index_from_ij( std::size_t const row, std::size_t const col ) const {
	std::size_t result;
	if( row < col ) {
		result = priv_index_from_ij( col, row );
	} else {
		result = row + col * ( 2 * m_num_vecs - col - 1 ) / 2;
	}
	return result;
}
}

