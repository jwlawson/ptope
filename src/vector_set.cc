/*
 * vector_set.cc
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
#include "vector_set.h"

namespace ptope {
template<class eT>
VectorSet<eT>::VectorSet( uint32_t const dimension, arma::uword initial_cap )
	: m_dimension { dimension }
	, m_ordered_pointers { VecPtrComparator( dimension ) }
	, m_vector_store ( dimension , initial_cap )
	, m_current_data_ptr { m_vector_store.memptr() }
{}
template<class eT>
bool
VectorSet<eT>::add( elem_t const * vec_ptr ) {
	std::size_t cur_size = m_ordered_pointers.size();
	if( cur_size == m_vector_store.n_cols ) {
		// Need to extend vector store
		priv_resize_extend();
	}
	return priv_insert_without_resize( vec_ptr );
}
template<class eT>
bool
VectorSet<eT>::priv_insert_without_resize( elem_t const * vec_ptr ) {

	bool inserted = false;
	auto eq_pair = m_ordered_pointers.equal_range(vec_ptr);

	// If the upper and lower bounds are different then the vector was found,
	// otherwise it is not in the store
	if(eq_pair.first == eq_pair.second) {
		auto& insert_hint = eq_pair.first;
		std::size_t cur_size = m_ordered_pointers.size();
		elem_t * next_store_vec = ptr_at( cur_size );
		std::memcpy( next_store_vec , vec_ptr , m_dimension * sizeof(elem_t) );
		m_ordered_pointers.insert( insert_hint , next_store_vec );
		inserted = true;
	}
	return inserted;
}
template<class eT>
void
VectorSet<eT>::priv_resize_extend() {
	arma::uword const new_size = m_vector_store.n_cols * 3 / 2 + 100;
	m_vector_store.resize( m_dimension , new_size );
	elem_t const * old_ptr = m_current_data_ptr;
	elem_t * new_ptr = m_vector_store.memptr();
	std::for_each( m_ordered_pointers.begin() , m_ordered_pointers.end() , 
			[old_ptr, new_ptr](elem_t const * const& p) {
				elem_t const*& q = const_cast<elem_t const* &>(p);
				q = new_ptr + std::distance( old_ptr , p );
			} );
	m_current_data_ptr = new_ptr;
}

template class VectorSet<double>;
template class VectorSet<float>;
template class VectorSet<uint32_t>;
template class VectorSet<uint64_t>;
template class VectorSet<int32_t>;
template class VectorSet<int64_t>;

}

