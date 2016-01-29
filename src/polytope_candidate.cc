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

#include "calc.h"

namespace ptope {
namespace {
constexpr double error = 10e-10;
/* As the program is never run in parallel with shared resources, we can safely
 * cache these at a program/global level */
arma::vec __new_vec_cached;
arma::vec __null_vec_cached;
arma::mat __qr_matrix;
arma::mat __solve_a;
arma::podarray<double> __solve_work;
arma::podarray<double> __qr_tau;
arma::podarray<double> __qr_work;
/* Simplification of arma::solve */
bool
solve(const arma::mat & A, const arma::vec & B) {
	using arma::uword;
	using arma::blas_int;
	__solve_a = A;
	const uword A_n_rows = A.n_rows;
	const uword A_n_cols = A.n_cols;

	const uword B_n_rows = B.n_rows;
	const uword B_n_cols = B.n_cols;

	char trans = 'N';
	blas_int  m     = blas_int(A_n_rows);
	blas_int  n     = blas_int(A_n_cols);
	blas_int  lda   = blas_int(A_n_rows);
	blas_int  ldb   = blas_int(A_n_cols);
	blas_int  nrhs  = blas_int(B_n_cols);
	blas_int  lwork = 3 * ( (std::max)(blas_int(1), m + (std::max)(m,nrhs)) );
	blas_int  info  = 0;

	__new_vec_cached.set_size(A_n_cols, B_n_cols);

	for(uword col=0; col<B_n_cols; ++col) {
		double* tmp_colmem = __new_vec_cached.colptr(col);
		arma::arrayops::copy( tmp_colmem, B.colptr(col), B_n_rows );
		for(uword row=B_n_rows; row<A_n_cols; ++row) {
			tmp_colmem[row] = double(0);
		}
	}
	__solve_work.set_min_size(static_cast<uword>(lwork));

	arma::lapack::gels<double>( &trans, &m, &n, &nrhs, __solve_a.memptr(), &lda,
			__new_vec_cached.memptr(), &ldb, __solve_work.memptr(), &lwork, &info );

	return (info == 0);
}
/** Simplification of arma::qr. We don't actually care about the result for R,
 * we only want to find the nullspace, so we can throw away a bunch of stuff. */
bool
qr_on_trans(const arma::mat & X) {
	using arma::uword;
	using arma::blas_int;
	/* This makes QR too big for just R, but the right size for Q. qeqrf fills QR
	 * with R in upper triangle and Q encoded in the bottom. orgqr decodes Q into
	 * the fill Q matrix. By Using the same matrix for everything we avoid having
	 * to copy all the values between Q and R between the function calls. */
	__qr_matrix.set_size(X.n_cols, X.n_cols);
	__qr_matrix.head_cols(X.n_rows) = X.t();

	const uword R_n_rows = __qr_matrix.n_rows;
	const uword R_n_cols = __qr_matrix.n_cols - 1;

	blas_int m         = static_cast<blas_int>(R_n_rows);
	blas_int n         = static_cast<blas_int>(R_n_cols);
	blas_int lwork     = 0;
	// take into account requirements of geqrf() _and_ orgqr()/ungqr()
	blas_int lwork_min = (std::max)(blas_int(1), (std::max)(m,n));
	blas_int k         = (std::min)(m,n);
	blas_int info      = 0;

	__qr_tau.set_min_size( static_cast<uword>(k) );

	double work_query[2];
	blas_int lwork_query = -1;

	arma::lapack::geqrf(&m, &n, __qr_matrix.memptr(), &m, __qr_tau.memptr(),
			&work_query[0], &lwork_query, &info);

	if(info != 0)  { return false; }

	blas_int lwork_proposed = static_cast<blas_int>(
			arma::access::tmp_real(work_query[0]) );
	lwork = (std::max)(lwork_proposed, lwork_min);
	__qr_work.set_min_size( static_cast<uword>(lwork) );

	arma::lapack::geqrf(&m, &n, __qr_matrix.memptr(), &m, __qr_tau.memptr(),
			__qr_work.memptr(), &lwork, &info);

	if(info != 0)  { return false; }

	arma::lapack::orgqr(&m, &m, &k, __qr_matrix.memptr(), &m, __qr_tau.memptr(),
			__qr_work.memptr(), &lwork, &info);
	return (info == 0);
}
}
PolytopeCandidate PolytopeCandidate::InValid;
PolytopeCandidate::PolytopeCandidate()
: _gram(std::make_shared<arma::mat>()),
	_vectors(arma::mat()),
	_basis_vecs_trans(),
	_hyperbolic(false),
	_valid(false) {}

PolytopeCandidate::PolytopeCandidate(const PolytopeCandidate & p)
: _gram(std::make_shared<arma::mat>(*p._gram)),
	_vectors(p._vectors),
	_basis_vecs_trans(p._basis_vecs_trans),
	_hyperbolic(p._hyperbolic),
	_valid(p._valid) {}

PolytopeCandidate::PolytopeCandidate(const GramMatrix & matrix)
: _gram(matrix),
	_vectors(arma::chol(*matrix)),
	_basis_vecs_trans(_vectors.underlying_matrix().t()),
	_hyperbolic(false),
	_valid(true) {}

PolytopeCandidate::PolytopeCandidate(GramMatrix && matrix)
: _gram(matrix),
	_vectors(arma::chol(*_gram)),
	_basis_vecs_trans(_vectors.underlying_matrix().t()),
	_hyperbolic(false),
	_valid(true) {}

PolytopeCandidate::PolytopeCandidate(const arma::mat & matrix)
: _gram(std::make_shared<arma::mat>(matrix)),
	_vectors(arma::chol(matrix)),
	_basis_vecs_trans(_vectors.underlying_matrix().t()),
	_hyperbolic(false),
	_valid(true) {}

PolytopeCandidate::PolytopeCandidate(arma::mat && matrix)
: _gram(std::make_shared<arma::mat>(matrix)),
	_vectors(arma::chol(*_gram)),
	_basis_vecs_trans(_vectors.underlying_matrix().t()),
	_hyperbolic(false),
	_valid(true) {}

PolytopeCandidate::PolytopeCandidate(const double * gram_ptr, int gram_size,
		const double * vector_ptr, int vector_dim, int no_vectors)
: _gram(std::make_shared<arma::mat>(gram_ptr, gram_size, gram_size)),
	_vectors(vector_ptr, vector_dim, no_vectors),
	_basis_vecs_trans(_vectors.underlying_matrix().t()),
	_hyperbolic(true),
	_valid(true) {}

PolytopeCandidate::PolytopeCandidate(
		std::initializer_list<std::initializer_list<double>> l)
: _gram(std::make_shared<arma::mat>(l)),
	_vectors(arma::chol(*_gram)),
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
PolytopeCandidate::vector_from_inner_products(const arma::vec & inner_vector)
		const {
	if(_hyperbolic) {
		solve(_basis_vecs_trans, inner_vector);
	} else {
		arma::solve(__new_vec_cached, _basis_vecs_trans, inner_vector);
	}
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
		qr_on_trans(_basis_vecs_trans);
		__null_vec_cached = __qr_matrix.col(__qr_matrix.n_cols - 1);
		const double ax = calc::mink_inner_prod(__null_vec_cached, __new_vec_cached);
		const double aa = calc::mink_sq_norm(__null_vec_cached);
		const double xx = calc::mink_sq_norm(__new_vec_cached);
		const double l = (-ax + std::sqrt(ax * ax + aa * (1 - xx)) )/aa;
		__null_vec_cached *= l;
		__new_vec_cached += __null_vec_cached;
	} else {
		const double e_norm = calc::eucl_sq_norm(__new_vec_cached);
		if(e_norm - 1.0 < error) {
			/* Invalid set of angles. */
			return false;
		}
		if(_hyperbolic) {
			__new_vec_cached(__new_vec_cached.size()-1) = std::sqrt(e_norm - 1.0);
		} else {
			const arma::uword last_entry = __new_vec_cached.size();
			__new_vec_cached.insert_rows(last_entry, 1, false);
			__new_vec_cached(last_entry) = std::sqrt(e_norm - 1.0);
		}
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
	result._gram->set_size(_gram->n_rows + 1, _gram->n_cols + 1);
	result._gram->submat(0, 0, _gram->n_rows - 1, _gram->n_cols - 1) = *_gram;
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

	const arma::uword last_col = _gram->n_cols;
	const arma::uword last_row = _gram->n_rows;
	for(arma::uword i = 0, max = result._vectors.size() - 1; i < max; ++i) {
		const double * const old_vec_ptr =
			result._vectors.underlying_matrix().colptr(i);
		const double val = calc::mink_inner_prod(new_vec.size(), new_vec.memptr(),
				old_vec_ptr);
		result._gram->at(i, last_col) = val;
		result._gram->at(last_row, i) = val;
	}
	result._gram->at(last_row, last_col) = calc::mink_sq_norm(new_vec);
}
void
PolytopeCandidate::rebase_vectors(arma::uvec vec_indices) {
	std::sort(vec_indices.begin(), vec_indices.end());
	for(arma::uword i = 0; i < vec_indices.size(); ++i) {
		_vectors.swap(i, vec_indices(i));
		_gram->swap_cols(i, vec_indices(i));
		_gram->swap_rows(i, vec_indices(i));
	}
	_basis_vecs_trans = _vectors.first_basis_cols().t();
	_basis_vecs_trans.unsafe_col(real_dimension()) *= -1;
}
PolytopeCandidate
PolytopeCandidate::swap_rebase(const arma::uword & a,
		const arma::uword & b) const {
	PolytopeCandidate result(*this);
	result._gram->swap_cols(a, b);
	result._gram->swap_rows(a, b);
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
	arma::vec evalues = arma::eig_sym(*_gram);
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
	_gram->save(os, arma::file_type::arma_binary);
	_vectors.save(os, arma::file_type::arma_binary);
	os << _hyperbolic << " ";
	os << _valid;
}
void
PolytopeCandidate::load(std::istream & is) {
	_gram->load(is, arma::file_type::arma_binary);
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
}
PolytopeCandidate &
PolytopeCandidate::operator=(const PolytopeCandidate & p) {
	*_gram = *p._gram;
	_basis_vecs_trans = p._basis_vecs_trans;
	_vectors = p._vectors;
	_hyperbolic = p._hyperbolic;
	_valid = p._valid;
	return *this;
}
std::ostream &
operator<<(std::ostream & os, const PolytopeCandidate & poly) {
	poly._gram->print(os, "Gram:");
	poly._vectors.underlying_matrix().print(os, "Vectors:");
	return os;
}
}

