template<class eT>
template<class T>
class VectorSet<eT>::vector_iterator {
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
	vector_iterator(elem_t * ptr, difference_type dim)
		: m_ptr(ptr), m_dimension(dim) {};

	reference operator*() const {
		m_returned = vec_t( m_ptr, m_dimension , false , true );
		return m_returned;
	}
	pointer   operator->() const
	{ return &( operator*() ); }
	reference operator[]( difference_type offset ) const {
		m_returned = vec_t( m_ptr + offset * m_dimension , m_dimension , false , true );
		return m_returned;
	}

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
	elem_t * m_ptr;
	difference_type m_dimension;
	mutable vec_t m_returned;
};

template<class eT>
VectorSet<eT>::VecPtrComparator::VecPtrComparator( uint32_t const size )
	: m_size(size)
{}
template<>
inline
bool
VectorSet<double>::VecPtrComparator::operator()( double const * lhs,
		double const * rhs ) const {
	static ptope::comparator::DoubleLess double_comp;
	return std::lexicographical_compare( lhs , lhs + m_size , rhs ,
			rhs + m_size , double_comp );
}
template<class eT>
inline
bool
VectorSet<eT>::VecPtrComparator::operator()( elem_t const * lhs,
		elem_t const * rhs ) const {
	return std::lexicographical_compare( lhs , lhs + m_size , rhs ,
			rhs + m_size );
}

template<class eT>
VectorSet<eT>::VectorSet( uint32_t const dimension, arma::uword initial_cap )
	: m_dimension { dimension }
	, m_ordered_pointers { VecPtrComparator( dimension ) }
	, m_vector_store ( dimension , initial_cap )
	, m_current_data_ptr { m_vector_store.memptr() }
{}
template<class eT>
inline
void
VectorSet<eT>::clear() {
	m_ordered_pointers.clear();
}
template<class eT>
inline
std::size_t
VectorSet<eT>::size() const {
	return m_ordered_pointers.size();
}
template<class eT>
inline
uint32_t
VectorSet<eT>::dimension() const {
	return m_vector_store.n_rows;
}
template<class eT>
inline
typename VectorSet<eT>::elem_t *
VectorSet<eT>::memptr() {
	return m_current_data_ptr;
}
template<class eT>
inline
typename VectorSet<eT>::elem_t const *
VectorSet<eT>::memptr() const {
	return m_current_data_ptr;
}
template<class eT>
inline
bool
VectorSet<eT>::contains( elem_t const * vec_ptr ) const {
	return m_ordered_pointers.find( vec_ptr ) != m_ordered_pointers.end();
}
template<class eT>
inline
bool
VectorSet<eT>::contains( vec_t const& vec ) const {
	return contains( vec.memptr() );
}
template<class eT>
inline
bool
VectorSet<eT>::add( vec_t const& vec) {
	return add( vec.memptr() );
}
template<class eT>
inline
typename VectorSet<eT>::vec_t
VectorSet<eT>::at( index_t const index ) {
	vec_t result( ptr_at(index) , m_dimension , false , true );
	return result;
}
template<class eT>
inline
typename VectorSet<eT>::elem_t *
VectorSet<eT>::ptr_at( index_t const index ) {
	return memptr() + index * m_dimension;
}
template<class eT>
inline
typename VectorSet<eT>::elem_t const *
VectorSet<eT>::ptr_at( index_t const index ) const {
	return memptr() + index * m_dimension;
}
template<class eT>
inline
typename VectorSet<eT>::iterator
VectorSet<eT>::begin()
{ return iterator( m_current_data_ptr , m_dimension ); }
template<class eT>
inline
typename VectorSet<eT>::const_iterator
VectorSet<eT>::begin() const
{ return const_iterator( m_current_data_ptr , m_dimension ); }
template<class eT>
inline
typename VectorSet<eT>::iterator
VectorSet<eT>::end()
{ return iterator( m_current_data_ptr + m_dimension * size(), m_dimension ); }
template<class eT>
inline
typename VectorSet<eT>::const_iterator
VectorSet<eT>::end() const
{ return const_iterator( m_current_data_ptr + m_dimension * size(), m_dimension ); }
