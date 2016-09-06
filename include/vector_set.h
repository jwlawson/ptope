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
//#include <boost/container/flat_set.hpp>

#include "armadillo"

namespace ptope {
class VectorSet {
	/** Compare the vectors which the pointers point to, using lexicographic
	 * ordering. */
	struct VecPtrComparator {
		VecPtrComparator( arma::uword const& size );
		bool
		operator()( double const * const lhs, double const * const rhs ) const;
	private:
		arma::uword const m_size;
	};
	typedef std::set<double const *, VecPtrComparator, boost::fast_pool_allocator<double const *>> PointerSet;
//	typedef boost::container::flat_set<double const *, VecPtrComparator> PointerSet;
	/** Iterator over the vectors. */
	template<class T>
	class vector_iterator {
	public:
		typedef vector_iterator<T> self_type;
		typedef std::ptrdiff_t difference_type;
		typedef T value_type;
		typedef T* pointer;
		typedef T& reference;
		typedef std::random_access_iterator_tag iterator_category;

		vector_iterator() = default;
		vector_iterator(self_type const&) = default;
		vector_iterator(self_type &&) = default;
		vector_iterator(double * ptr, difference_type dim);

		reference operator*() const;
		pointer   operator->() const;
		reference operator[]( difference_type offset ) const;

		self_type & operator++()
		{ m_ptr += m_dimension; return *this; }
		self_type   operator++( int )
		{ self_type tmp(*this); m_ptr += m_dimension; return tmp; }
		self_type & operator--()
		{ m_ptr -= m_dimension; return *this; }
		self_type   operator--( int )
		{ self_type tmp(*this); m_ptr -= m_dimension; return tmp; }

		self_type & operator+=( difference_type offset )
		{ m_ptr += offset * m_dimension; return *this; }
		self_type & operator-=( difference_type offset )
		{ m_ptr -= offset * m_dimension; return *this; }

		friend self_type operator+( self_type other , difference_type offset )
		{ other += offset; return other; }
		friend self_type operator+( difference_type offset , self_type other )
		{ other += offset; return other; }
		friend self_type operator-( self_type other , difference_type offset )
		{ other -= offset; return other; }
		friend self_type operator-( difference_type offset , self_type other )
		{ other -= offset; return other; }
		friend difference_type operator-( self_type const& lhs , self_type const& rhs )
		{ return ( lhs.m_ptr - rhs.m_ptr ) / lhs.m_size; }

		friend bool operator==( self_type const& lhs , self_type const& rhs )
		{ return lhs.m_ptr == rhs.m_ptr; }
		friend bool operator!=( self_type const& lhs , self_type const& rhs )
		{ return !operator==(lhs, rhs); }
		friend bool operator< ( self_type const& lhs , self_type const& rhs )
		{ return rhs - lhs > 0; }
		friend bool operator> ( self_type const& lhs , self_type const& rhs )
		{ return rhs < lhs; }
		friend bool operator<=( self_type const& lhs , self_type const& rhs )
		{ return !(lhs > rhs); }
		friend bool operator>=( self_type const& lhs , self_type const& rhs )
		{ return !(lhs < rhs); }
	private:
		double * m_ptr;
		difference_type m_dimension;
		mutable arma::vec m_returned;
	};
public:
	typedef vector_iterator<arma::vec> iterator;
	typedef vector_iterator<arma::vec const> const_iterator;

	VectorSet( arma::uword const& dimension, arma::uword initial_cap = 10 );

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
	std::size_t dimension() const;
	/**
	 * Check whether the provided vector is already in the set.
	 * Return: true if present in set
	 * Complexity: logarithmic
	 */
	bool contains( double const * vec_ptr ) const;
	bool contains( arma::vec const& vec ) const;
	/**
	 * Insert the provided vector into the set.
	 * Return: true if inserted, false if already present
	 * Complexity: Logarithmic search plus linear insertion
	 */
	bool add( double const * vec_ptr );
	bool add( arma::vec const& vec);
	/**
	 * Get a (strict) vector from the set at the specified index.
	 * This vector uses the memory in the set, so cannot be resized and should not
	 * be changed.
	 *
	 * Also, if the vector store is resized, this vector will be invalidated.
	 */
	arma::vec at( arma::uword const& index );
	double * ptr_at( arma::uword const& index );
	double const * ptr_at( arma::uword const& index ) const;
	/**
	 * Get a pointer to the underlying matrix of vectors.
	 */
	double * memptr();
	double const * memptr() const;
	/**
	 * Get iterators to the vectors.
	 */
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;

private:
	arma::uword const m_dimension;
	PointerSet m_ordered_pointers;
	arma::mat m_vector_store;
	double * m_current_data_ptr;

	bool
	priv_insert_without_resize( double const * vec_ptr );
	void
	priv_resize_extend();
};


inline
void
VectorSet::clear() {
	m_ordered_pointers.clear();
}
inline
std::size_t
VectorSet::size() const {
	return m_ordered_pointers.size();
}
inline
std::size_t
VectorSet::dimension() const {
	return m_vector_store.n_rows;
}
inline
double *
VectorSet::memptr() {
	return m_current_data_ptr;
}
inline
double const *
VectorSet::memptr() const {
	return m_current_data_ptr;
}
inline
arma::vec
VectorSet::at( arma::uword const& index ) {
	arma::vec result( ptr_at(index) , m_dimension , false , true );
	return result;
}
inline
double *
VectorSet::ptr_at( arma::uword const& index ) {
	return memptr() + index * m_dimension;
}
inline
double const *
VectorSet::ptr_at( arma::uword const& index ) const {
	return memptr() + index * m_dimension;
}
inline
VectorSet::iterator
VectorSet::begin()
{ return iterator( m_current_data_ptr , m_dimension ); }
inline
VectorSet::const_iterator
VectorSet::begin() const
{ return const_iterator( m_current_data_ptr , m_dimension ); }
inline
VectorSet::iterator
VectorSet::end()
{ return iterator( m_current_data_ptr + m_dimension * size(), m_dimension ); }
inline
VectorSet::const_iterator
VectorSet::end() const
{ return const_iterator( m_current_data_ptr + m_dimension * size(), m_dimension ); }

template<class T>
VectorSet::vector_iterator<T>::vector_iterator(double * ptr, difference_type dim)
	: m_ptr(ptr), m_dimension(dim) {}

template<class T>
typename VectorSet::vector_iterator<T>::reference
VectorSet::vector_iterator<T>::operator*() const {
	m_returned = arma::vec( m_ptr, m_dimension , false , true );
	return m_returned;
}
template<class T>
typename VectorSet::vector_iterator<T>::pointer
VectorSet::vector_iterator<T>::operator->() const
{ return &( operator*() ); }

template<class T>
typename VectorSet::vector_iterator<T>::reference
VectorSet::vector_iterator<T>::operator[]( difference_type offset ) const {
	m_returned = arma::vec( m_ptr + offset * m_dimension , m_dimension , false , true );
	return m_returned;
}
}
#endif

