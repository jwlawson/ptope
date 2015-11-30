/*
 * polytope_candidate.h
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
#include "polytope_candidate.h"

namespace ptope {
namespace {
constexpr double error = 1e-14;
double mink_inner_prod(const arma::vec & a, const arma::vec & b) {
	double sq = 0;
	arma::uword max = a.size() - 1;
	for(arma::uword i = 0; i < max; ++i) {
		sq += a(i) * b(i);
	}
	sq -= a(max) * b(max);
	return sq;
}
double mink_inner_prod(const arma::vec & a) {
	double sq = 0;
	arma::uword max = a.size() - 1;
	for(arma::uword i = 0; i < max; ++i) {
		sq += a(i) * a(i);
	}
	sq -= a(max) * a(max);
	return sq;
}
double eucl_sq_norm(const arma::vec & a) {
	double sq = 0;
	arma::uword max = a.size();
	for(arma::uword i = 0; i < max; ++i) {
		sq += a(i) * a(i);
	}
	return sq;
}
}
PolytopeCandidate PolytopeCandidate::InValid;
PolytopeCandidate::PolytopeCandidate(const GramMatrix & matrix)
	: _gram(matrix),
		_vectors(arma::chol(matrix)),
		_basis_vecs_trans(_vectors.underlying_matrix().t()),
		_hyperbolic(false),
		_valid(true) {}

PolytopeCandidate::PolytopeCandidate(GramMatrix && matrix)
	: _gram(matrix),
		_vectors(arma::chol(matrix)),
		_basis_vecs_trans(_vectors.underlying_matrix().t()),
		_hyperbolic(false),
		_valid(true) {}

PolytopeCandidate
PolytopeCandidate::extend_by_inner_products(const arma::vec & inner_vector) const {
	arma::vec new_vec = arma::solve(_basis_vecs_trans, inner_vector);
	if(_hyperbolic && std::abs(new_vec(new_vec.size() - 1)) > error) {
		/* 
		 * Rescale the new vector by adding something from the nullspace, so that
		 * the norm of the vector is 1.
		 *
		 * The 'basis' vectors are not actually a basis, but leave one dimension
		 * undefined. This is the nullspace vector n here.
		 *
		 * As this is the nullspace, if x is any vector satisfying Ax = b, then
		 * (x + l*a) will also satisfy the same equation. Hence to fins a unit
		 * vector satisfying the equation we need to find l such that (x + l*a) is
		 * unit, or <x+la, x+la> = 1. This simplifies to a quadratic equation:
		 * 	1 = <x,x> + 2l<a,x> + l*l<a,a>
		 * which has solution:
		 * 	l = (-<a,x> + sqrt( <a,x>*<a,x> + <a,a>(<x,x> - 1) ))/<a,a>
		 */
		// arma::null was introduced in arma version 5
		arma::vec n = arma::null(_basis_vecs_trans);
		const double ax = mink_inner_prod(n, new_vec);
		const double aa = mink_inner_prod(n);
		const double xx = mink_inner_prod(new_vec);
		const double l = (-ax + std::sqrt(ax * ax + aa * (1 - xx)) )/aa;
		new_vec = new_vec + (l * n);
	} else {
		const double e_norm = eucl_sq_norm(new_vec);
		if(e_norm - 1.0 < error) {
			/* Invalid set of angles. */
			return PolytopeCandidate::InValid;
		}
		if(_hyperbolic) {
			new_vec(new_vec.size()-1) = std::sqrt(e_norm - 1.0);
		} else {
			const arma::uword last_entry = new_vec.size();
			new_vec.insert_rows(last_entry, 1, false);
			new_vec(last_entry) = std::sqrt(e_norm - 1.0);
		}
	}
	PolytopeCandidate result(*this);
	/* The copy allocates memory for the matrix, then this resize call will also
	 * allocate and copy the memory. There is probably a better way to avoid this
	 * double allocation. */
	result._gram.resize(_gram.n_rows + 1, _gram.n_cols + 1);
	if(!_hyperbolic) {
		result._vectors.add_first_hyperbolic_vector(new_vec);
		result._basis_vecs_trans = result._vectors.first_basis_cols().t();
		result._hyperbolic = true;
	} else {
		result._vectors.add_vector(new_vec);
	}
	/* Add new inner products to the gram matrix */
	const arma::uword last_col = _gram.n_cols;
	const arma::uword last_row = _gram.n_rows;
	for(arma::uword i = 0, max = result._vectors.size() - 1;
			i < max; ++i) {
		result._gram(i, last_col) =
			mink_inner_prod(new_vec, result._vectors.unsafe_get(i));
		result._gram(last_row, i) =
			mink_inner_prod(new_vec, result._vectors.unsafe_get(i));
	}
	result._gram(last_row, last_col) = mink_inner_prod(new_vec);
	return result;
}
void
PolytopeCandidate::rebase_vectors(arma::uvec vec_indices) {
	std::sort(vec_indices.begin(), vec_indices.end());
	for(arma::uword i = 0; i < vec_indices.size(); ++i) {
		_vectors.swap(i, vec_indices(i));
		_gram.swap_cols(i, vec_indices(i));
		_gram.swap_rows(i, vec_indices(i));
	}
	_basis_vecs_trans = _vectors.first_basis_cols().t();
	_basis_vecs_trans.unsafe_col(real_dimension()) *= -1;
}
PolytopeCandidate
PolytopeCandidate::swap_rebase(const arma::uword & a,
		const arma::uword & b) const {
	PolytopeCandidate result(*this);
	result._gram.swap_cols(a, b);
	result._gram.swap_rows(a, b);
	result._vectors.swap(a, b);
	result._basis_vecs_trans = _vectors.first_basis_cols().t();
	result._basis_vecs_trans.unsafe_col(real_dimension()) *= -1;
	return result;
}
void
PolytopeCandidate::recompute_gram() {
	arma::uword cm = _gram.n_cols;
	for(arma::uword row = 0, rm = _gram.n_rows; row < rm; ++row) {
		for(arma::uword col = 0; col < cm; ++col) {
			_gram(row,col) = mink_inner_prod(_vectors.unsafe_get(row),
					_vectors.unsafe_get(col));
		}
	}
}
std::pair<uint,uint>
PolytopeCandidate::signature() const {
	arma::vec evalues = arma::eig_sym(_gram);
	auto result = std::make_pair((uint)0, (uint) 0);
	for(arma::uword i = 0; i < evalues.size(); ++i) {
		if (std::abs(evalues(i)) < error) {
			/* Zero determinant */
		} else if(evalues(i) < 0) {
			++(result.second);
		} else {
			++(result.first);
		}
	}
	return result;
}
std::size_t
PolytopeCandidate::real_dimension() const {
	if(_hyperbolic) {
		return static_cast<std::size_t>(_vectors.dimension() - 1);
	} else {
		return static_cast<std::size_t>(_vectors.dimension());
	}
}
void
PolytopeCandidate::save(std::ostream & os) {
	_gram.save(os, arma::file_type::arma_binary);
	_vectors.save(os, arma::file_type::arma_binary);
	os << _hyperbolic << " ";
	os << _valid;
}
void
PolytopeCandidate::load(std::istream & is) {
	_gram.load(is, arma::file_type::arma_binary);
	_vectors.load(is, arma::file_type::arma_binary);
	_basis_vecs_trans = _vectors.first_basis_cols().t();
	is >> _hyperbolic;
	is >> _valid;
}
std::ostream &
operator<<(std::ostream & os, const PolytopeCandidate & poly) {
	poly._gram.print(os, "Gram:");
	poly._vectors.underlying_matrix().print(os, "Vectors:");
	return os;
}
}

