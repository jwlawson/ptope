/*
 * polytope_candidate.cc
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
constexpr double error = 1e-15;
double mink_inner_prod(const arma::vec & a, const arma::vec & b) {
	double sq = 0;
	arma::uword max = a.size() - 1;
	for(arma::uword i = 0; i < max; ++i) {
		sq += a(i) * b(i);
	}
	sq -= a(max) * b(max);
	return std::sqrt(sq);
}
double mink_inner_prod(const arma::vec & a) {
	double sq = 0;
	arma::uword max = a.size() - 1;
	for(arma::uword i = 0; i < max; ++i) {
		sq += a(i) * a(i);
	}
	sq -= a(max) * a(max);
	return std::sqrt(sq);
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
PolytopeCandidate::extend_by_inner_products(const arma::vec & inner_vector) {
	arma::vec new_vec = arma::solve(_basis_vecs_trans, inner_vector);
	const double e_norm = eucl_sq_norm(new_vec);
	if(e_norm < 1) {
		/* Invalid set of angles. */
		return PolytopeCandidate::InValid;
	}
	PolytopeCandidate result(*this);
	/* The copy allocates memory for the matrix, then this resize call will also
	 * allocate and copy the memory. There is probably a better way to avoid this
	 * double allocation. */
	result._gram.resize(_gram.n_rows + 1, _gram.n_cols + 1);
	if(!_hyperbolic) {
		const arma::uword last_entry = new_vec.size();
		new_vec.insert_rows(last_entry, 1, false);
		new_vec(last_entry) = std::sqrt(e_norm - 1);
		result._vectors.add_first_hyperbolic_vector(new_vec);
		result._basis_vecs_trans = result._vectors.first_basis_cols().t();
		result._hyperbolic = true;
	} else {
		new_vec(new_vec.size()-1) = std::sqrt(e_norm - 1);
		result._vectors.add_vector(new_vec);
	}
	/* Add new inner products to the gram matrix */
	const arma::uword last_col = _gram.n_cols;
	const arma::uword last_row = _gram.n_rows;
	for(arma::uword i = 0; i < inner_vector.size(); ++i) {
		result._gram(i, last_col) = inner_vector(i);
		result._gram(last_row, i) = inner_vector(i);
	}
	for(arma::uword i = inner_vector.size(), max = result._vectors.size() - 1;
			i < max; ++i) {
		result._gram(i, last_col) = mink_inner_prod(new_vec, _vectors.unsafe_get(i));
		result._gram(last_row, i) = mink_inner_prod(new_vec, _vectors.unsafe_get(i));
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
std::ostream &
operator<<(std::ostream & os, const PolytopeCandidate & poly) {
	poly._gram.print(os, "Gram:");
	poly._vectors.underlying_matrix().print(os, "Vectors:");
	return os;
}
}

