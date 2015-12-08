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
/* As the program is never run in parallel with shared resources, we can safely
 * cache these at a program/global level */
arma::vec __new_vec_cached;
arma::vec __null_vec_cached;
#ifdef VALGRIND_SAFE
arma::mat __q_cached;
arma::mat __r_cached;
#endif
}
PolytopeCandidate PolytopeCandidate::InValid;
PolytopeCandidate::PolytopeCandidate()
: _gram(),
	_vectors(arma::mat()),
	_basis_vecs_trans(),
	_hyperbolic(false),
	_valid(false) {}

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

PolytopeCandidate::PolytopeCandidate(const double * gram_ptr, int gram_size,
		const double * vector_ptr, int vector_dim, int no_vectors)
: _gram(gram_ptr, gram_size, gram_size),
	_vectors(vector_ptr, vector_dim, no_vectors),
	_basis_vecs_trans(_vectors.underlying_matrix().t()),
	_hyperbolic(true),
	_valid(true) {}

PolytopeCandidate
PolytopeCandidate::extend_by_inner_products(const arma::vec & inner_vector) const {
	arma::solve(__new_vec_cached, _basis_vecs_trans, inner_vector);
	if(_hyperbolic && 
			std::abs(__new_vec_cached(__new_vec_cached.size() - 1)) > error) {
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
#ifdef VALGRIND_SAFE
		arma::qr(__q_cached, __r_cached, _basis_vecs_trans.t());
		__null_vec_cached = __q_cached.col(__q_cached.n_cols - 1);
#else
		// arma::null was introduced in arma version 5 and not supported by Valgrind
		// 3.10.1
		arma::null(__null_vec_cached, _basis_vecs_trans);
#endif
		const double ax = mink_inner_prod(__null_vec_cached, __new_vec_cached);
		const double aa = mink_inner_prod(__null_vec_cached);
		const double xx = mink_inner_prod(__new_vec_cached);
		const double l = (-ax + std::sqrt(ax * ax + aa * (1 - xx)) )/aa;
		__new_vec_cached += (l * __null_vec_cached);
	} else {
		const double e_norm = eucl_sq_norm(__new_vec_cached);
		if(e_norm - 1.0 < error) {
			/* Invalid set of angles. */
			return PolytopeCandidate::InValid;
		}
		if(_hyperbolic) {
			__new_vec_cached(__new_vec_cached.size()-1) = std::sqrt(e_norm - 1.0);
		} else {
			const arma::uword last_entry = __new_vec_cached.size();
			__new_vec_cached.insert_rows(last_entry, 1, false);
			__new_vec_cached(last_entry) = std::sqrt(e_norm - 1.0);
		}
	}
	return extend_by_vector(__new_vec_cached);
}
PolytopeCandidate
PolytopeCandidate::extend_by_vector(const arma::vec & new_vec) const {
	PolytopeCandidate result;
	result._gram.set_size(_gram.n_rows + 1, _gram.n_cols + 1);
	result._gram.submat(0, 0, _gram.n_rows - 1, _gram.n_cols - 1) = _gram;
	result._hyperbolic = true;
	result._valid = true;
	if(!_hyperbolic) {
		result._vectors.copy_and_add_first_hyperbolic_vector(_vectors, new_vec);
		result._basis_vecs_trans = result._vectors.first_basis_cols().t();
		/* Note, don't need to multiply last column by -1 as the values are all 0 */
	} else {
		result._vectors.copy_and_add_vector(_vectors, new_vec);
		result._basis_vecs_trans = _basis_vecs_trans;
	}

	const arma::uword last_col = _gram.n_cols;
	const arma::uword last_row = _gram.n_rows;
	for(arma::uword i = 0, max = result._vectors.size() - 1; i < max; ++i) {
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
	result._basis_vecs_trans = result._vectors.first_basis_cols().t();
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

