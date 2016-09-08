/*
 * vector_set.h
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
 * A set of vectors. As with other sets, this ensures that each vector in the
 * set is unique. The vectors are not sorted in any way, but are stored in a
 * consistent ordering given by how they were added.
 */
#pragma once
#ifndef _PTOPE_VECTOR_SET_H_
#define _PTOPE_VECTOR_SET_H_

#include <set>

#include <boost/pool/pool_alloc.hpp>

#include <armadillo>

#include "comparator.h"

namespace ptope {
template<class eT>
class VectorSet {
	/** Compare the vectors which the pointers point to, using lexicographic
	 * ordering. */
	struct VecPtrComparator {
		VecPtrComparator( uint32_t const size );
		bool operator()( eT const * const lhs, eT const * const rhs ) const;
	private:
		uint32_t const m_size;
	};

	template<class T>
	class vector_iterator;
public:
	typedef eT elem_t;
	typedef VectorSet<elem_t> self_t;
	typedef arma::Col<elem_t> vec_t;
	typedef arma::uword index_t;
	typedef vector_iterator<vec_t> iterator;
	typedef vector_iterator<vec_t const> const_iterator;
	typedef std::set<elem_t const *, VecPtrComparator,
					boost::fast_pool_allocator<elem_t const *>> PointerSet;

	VectorSet( uint32_t const dimension, arma::uword initial_cap = 10 );

	/**
	 * Clear the set.
	 *
	 * No memory is released.
	 */
	void clear();
	/**
	 * Get the number of vectors in the set.
	 */
	std::size_t size() const;
	/**
	 * Get the dimension of each of the vectors.
	 */
	uint32_t dimension() const;
	/**
	 * Check whether the provided vector is already in the set.
	 * Return: true if present in set
	 * Complexity: logarithmic
	 */
	bool contains( elem_t const * vec_ptr ) const;
	bool contains( vec_t const& vec ) const;
	/**
	 * Insert the provided vector into the set.
	 * Return: true if inserted, false if already present
	 * Complexity: Logarithmic search plus linear insertion
	 */
	bool add( elem_t const * vec_ptr );
	bool add( vec_t const& vec);
	/**
	 * Get a (strict) vector from the set at the specified index.
	 * This vector uses the memory in the set, so cannot be resized and should not
	 * be changed.
	 *
	 * Also, if the vector store is resized, this vector will be invalidated.
	 */
	vec_t at( index_t const index );
	elem_t * ptr_at( index_t const index );
	elem_t const * ptr_at( index_t const index ) const;
	/**
	 * Get a pointer to the underlying matrix of vectors.
	 */
	elem_t * memptr();
	elem_t const * memptr() const;
	/**
	 * Get iterators to the vectors.
	 */
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;

private:
	uint32_t const m_dimension;
	PointerSet m_ordered_pointers;
	arma::Mat<elem_t> m_vector_store;
	elem_t * m_current_data_ptr;

	bool
	priv_insert_without_resize( elem_t const * vec_ptr );
	void
	priv_resize_extend();
};

#include "detail/vector_set.inl"

}
#endif

