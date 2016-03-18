/*
 * lq_info.h
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
/**
 * Struct containing the LQ decomposition of a PolytopeCandidate's basis vector
 * transformation matrix. By keeping a copy of this decomposition it does not
 * need to be computed every time the PolytopeCandidate is extended.
 */
#pragma once
#ifndef PTOPE_LQ_INFO_H_
#define PTOPE_LQ_INFO_H_

#include <armadillo>
#include <memory>

namespace ptope {
namespace detail {
class LQInfo {
public:
	static
	std::unique_ptr<LQInfo>
	compute(arma::mat const & m);
	/** Get reference to Q**T * inv(L) matrix */
	arma::mat const &
	qtli() const;
	/** Get reference to nullspace vector */
	arma::vec const &
	null() const;
private:
	/**
	 * Construct the LQ decomposition of the provided matrix using LAPACK deglqf.
	 */
	LQInfo(arma::uword const size);
	/**
	 * Product Q**T * inv(L)
	 */
	arma::mat _qtli;
	/**
	 * Vector of nullspace of the system of equations defined by LQ.
	 */
	arma::vec _nullspace;
};
inline
arma::mat const &
LQInfo::qtli() const {
	return _qtli;
}
inline
arma::vec const &
LQInfo::null() const {
	return _nullspace;
}
#if !defined(ARMA_BLAS_CAPITALS)
#define arma_dorglq dorglq
#define arma_dgelqf dgelqf
#else
#define arma_dorglq DORGLQ
#define arma_dgelqf DGEQLF
#endif
extern "C" {
/* Compute LQ decomposition of matrix. */
void arma_fortran(arma_dgelqf)(
	arma::blas_int* m,
	arma::blas_int* n,
	double* a,
	arma::blas_int* lda,
	double* tau,
	double* work,
	arma::blas_int*	lwork,
	arma::blas_int* info);
/* Find Q from the LQ decom given by dgelqf */
void arma_fortran(arma_dorglq)(
	arma::blas_int* m,
	arma::blas_int* n,
	arma::blas_int* k,
	double* a,
	arma::blas_int* lda,
	double* tau,
	double* work,
	arma::blas_int* lwork,
	arma::blas_int* info);
}
}
}
#endif

