inline
void
GramMatrix::from( VectorSet<double> const& vectors ) {
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
double
GramMatrix::raw_rfp( std::size_t const& index ) const {
	return m_matrix[ index ];
}
inline
std::size_t
GramMatrix::priv_min_product_size( VectorSet<double> const& vectors ) {
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
