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

#include "calc.h"

namespace ptope {
namespace {
constexpr double error = 10e-10;
/* As the program is never run in parallel with shared resources, we can safely
 * cache these at a program/global level */
arma::vec __new_vec_cached;
arma::vec __null_vec_cached;
}
/* Static private vars */
PolytopeCandidate PolytopeCandidate::InValid;

PolytopeCandidate::PolytopeCandidate()
: _gram(),
	_vectors(arma::mat()),
	_basis_vecs_trans(),
	_hyperbolic(false),
	_valid(false) {}

PolytopeCandidate::PolytopeCandidate(const PolytopeCandidate & p)
: _gram(p._gram),
	_vectors(p._vectors),
	_basis_vecs_trans(p._basis_vecs_trans),
	_hyperbolic(p._hyperbolic),
	_valid(p._valid),
	_lq_info(p._lq_info ? new detail::LQInfo(*p._lq_info) : nullptr) {}

PolytopeCandidate::PolytopeCandidate(const arma::mat & matrix)
: _gram(matrix),
	_vectors(arma::chol(matrix)),
	_basis_vecs_trans(_vectors.underlying_matrix().t()),
	_hyperbolic(false),
	_valid(true) {}

PolytopeCandidate::PolytopeCandidate(arma::mat && matrix)
: _gram(std::move(matrix)),
	_vectors(arma::chol(_gram)),
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

PolytopeCandidate::PolytopeCandidate(
		std::initializer_list<std::initializer_list<double>> l)
: _gram(l),
	_vectors(arma::chol(_gram)),
	_basis_vecs_trans(_vectors.underlying_matrix().t()),
	_hyperbolic(false),
	_valid(true) {}

PolytopeCandidate
PolytopeCandidate::extend_by_inner_products(const arma::vec & inner_vector) const {
	if(vector_from_inner_products(inner_vector)) {
		return extend_by_vector(__new_vec_cached);
	} else {
		return PolytopeCandidate::InValid;
	}
}
bool
PolytopeCandidate::extend_by_inner_products(PolytopeCandidate & result,
		const arma::vec & inner_vector) const {
	if(vector_from_inner_products(inner_vector)) {
		extend_by_vector(result, __new_vec_cached);
		return true;
	} else {
		return false;
	}
}
bool
PolytopeCandidate::vector_from_inner_products(const arma::vec & inner_vector) const {
	if(_hyperbolic) {
		if(!_lq_info) {
			_lq_info = detail::LQInfo::compute(_basis_vecs_trans);
		}
		__null_vec_cached = _lq_info->null();
		__new_vec_cached = _lq_info->qtli() * inner_vector;
		/* 
		 * Rescale the new vector by adding something from the nullspace, so that
		 * the norm of the vector is 1.
		 *
		 * The 'basis' vectors are not actually a basis, but leave one dimension
		 * undefined. This is the nullspace vector n here.
		 *
		 * As this is the nullspace, if x is any vector satisfying Ax = b, then
		 * (x + l*a) will also satisfy the same equation. Hence to find a unit
		 * vector satisfying the equation we need to find l such that (x + l*a) is
		 * unit, or <x+la, x+la> = 1. This simplifies to a quadratic equation:
		 * 	1 = <x,x> + 2l<a,x> + l*l<a,a>
		 * which has solution:
		 * 	l = (-<a,x> + sqrt( <a,x>*<a,x> + <a,a>(<x,x> - 1) ))/<a,a>
		 */
		const double xx = calc::mink_sq_norm(__new_vec_cached);
		if(std::abs(xx - 1.0) < error) {
			/* Lapack solver returns a (Euclidean) unit vector. If the Minkowski norm
			 * is also 1, then we have constructed a non-hyperbolic vector, which we
			 * don't want. */
			return false;
		}
		const double ax = calc::mink_inner_prod(__null_vec_cached, __new_vec_cached);
		const double aa = calc::mink_sq_norm(__null_vec_cached);
		const double disc = ax * ax + aa * (1.0 - xx);
		if(disc < 0) {
			/* If discriminant is negative then no solutions */
			return false;
		}
		double l;
		if(std::abs(aa + 1) < error) {
			/* Nullspace is vector (0, 0, 0, ..., 1).
			 * The below also works, but we know that in case we want lm and that
			 * ax = 0, so the computation can be simplified quite a bit. */
			l = std::sqrt(disc);
		} else {
			double lm = (-ax - std::sqrt(disc) )/aa;
			double lp = (-ax + std::sqrt(disc) )/aa;
			bool plus = true;
			bool minus = true;
			for(arma::uword vec_ind = inner_vector.size(), max = _vectors.size();
					vec_ind < max && (plus || minus);
					++vec_ind) {
				const double * v = _vectors.get_ptr(vec_ind);
				const double xv = calc::mink_inner_prod(__new_vec_cached.size(),
						__new_vec_cached.memptr(), v);
				const double av = calc::mink_inner_prod(__new_vec_cached.size(),
						__null_vec_cached.memptr(), v);
				if(xv + lp * av > error) {
					plus = false;
				}
				if(xv + lm * av > error) {
					minus = false;
				}
			}
			if(plus) {
				l = lp;
			} else if(minus) {
				l = lm;
			} else {
				return false;
			}
		}
		__null_vec_cached *= l;
		__new_vec_cached += __null_vec_cached;
	} else {
		arma::solve(__new_vec_cached, _basis_vecs_trans, inner_vector);
		const double e_norm = calc::eucl_sq_norm(__new_vec_cached);
		if(e_norm - 1.0 < error) {
			/* Invalid set of angles. */
			return false;
		}
		const arma::uword last_entry = __new_vec_cached.size();
		__new_vec_cached.insert_rows(last_entry, 1, false);
		__new_vec_cached(last_entry) = std::sqrt(e_norm - 1.0);
	}
	return true;
}
PolytopeCandidate
PolytopeCandidate::extend_by_vector(const arma::vec & new_vec) const {
	PolytopeCandidate result;
	extend_by_vector(result, new_vec);
	return result;
}
void
PolytopeCandidate::extend_by_vector(PolytopeCandidate & result,
		const arma::vec & new_vec) const {
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
		const double * const old_vec_ptr = result._vectors.get_ptr(i);
		const double val = calc::mink_inner_prod(new_vec.size(), new_vec.memptr(),
				old_vec_ptr);
		result._gram.at(i, last_col) = val;
		result._gram.at(last_row, i) = val;
	}
	result._gram.at(last_row, last_col) = calc::mink_sq_norm(new_vec);
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
/* Because we force the vectors to live in the space with signature (d,1), the
 * matrix should always have this signature too. This means this calculation is
 * generally pointless - but could be useful to check that everything is working
 * as it should. */
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
PolytopeCandidate::save(std::ostream & os) const {
	_gram.save(os, arma::file_type::arma_binary);
	_vectors.save(os, arma::file_type::arma_binary);
	os << _hyperbolic << " ";
	os << _valid << '\n';
}
void
PolytopeCandidate::load(std::istream & is) {
	_gram.load(is, arma::file_type::arma_binary);
	_vectors.load(is, arma::file_type::arma_binary);
	_basis_vecs_trans = _vectors.first_basis_cols().t();
	is >> _hyperbolic;
	is >> _valid;
}
void
PolytopeCandidate::swap(PolytopeCandidate & p) {
	_gram.swap(p._gram);
	_vectors.swap(p._vectors);
	_basis_vecs_trans.swap(p._basis_vecs_trans);
	std::swap(_hyperbolic, p._hyperbolic);
	std::swap(_valid, p._valid);
	std::swap(_lq_info, p._lq_info);
}
std::ostream &
operator<<(std::ostream & os, const PolytopeCandidate & poly) {
	poly._gram.print(os, "Gram:");
	poly._vectors.underlying_matrix().print(os, "Vectors:");
	return os;
}
PolytopeCandidate &
PolytopeCandidate::operator=(const PolytopeCandidate & p) {
	_gram = p._gram;
	_vectors = p._vectors;
	_basis_vecs_trans = p._basis_vecs_trans;
	_hyperbolic = p._hyperbolic;
	_valid = p._valid;
	/* TODO Check if this is best way of handling assignment of LQ.
	 * Depends when LQ info is computed and when it's needed. */
	_lq_info.release();
	return *this;
}
}

