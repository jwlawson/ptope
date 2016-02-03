/*
 * underdetermined_solver.h
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
#pragma once
#ifndef _PTOPE_UNDERDETERMINED_SOLVER_H_
#define _PTOPE_UNDERDETERMINED_SOLVER_H_

#include <armadillo>

namespace ptope {
/** 
 * Solve an underdetermined system of linear equations, using BLAS/LAPACK dgels.
 */
struct UDSolver {
public:
	bool
	operator()(arma::vec & out, arma::vec & nullvec,
			const arma::mat & A, const arma::vec & B);
private:
	arma::mat _a;
	arma::podarray<double> _work;
	/** Cached tau work array for finding LQ.
	 * TODO Could be incorporated into _work, so that the two work arrays are
	 * always next to each other (hence cached properly) */
	arma::podarray<double> _tau;
};
}

#if !defined(ARMA_BLAS_CAPITALS)
#define arma_dorglq dorglq
#define arma_dgelqf dgelqf
#define arma_dormlq dormlq
#else
#define arma_dorglq DORGLQ
#define arma_dgelqf DGEQLF
#define arma_dormlq DORMLQ
#endif
extern "C" {
/* Compute LQ decomposition of matrix. */
	void arma_fortran(dgelqf)(
		arma::blas_int* m,
		arma::blas_int* n,
		double* a,
		arma::blas_int* lda,
		double* tau,
		double* work,
		arma::blas_int*	lwork,
		arma::blas_int* info);
/* Find Q from the LQ decom given by dgelqf */
	void arma_fortran(dorglq)(
		arma::blas_int* m,
		arma::blas_int* n,
		arma::blas_int* k,
		double* a,
		arma::blas_int* lda,
		double* tau,
		double* work,
		arma::blas_int* lwork,
		arma::blas_int* info);
	void arma_fortran(arma_dormlq)(
		char* side,
		char* trans,
		arma::blas_int* m,
		arma::blas_int* n,
		arma::blas_int* k,
		double* a,
		arma::blas_int* lda,
		double* tau,
		double* c,
		arma::blas_int* ldc,
		double* work,
		arma::blas_int* lwork,
		arma::blas_int* info);
}
#endif

