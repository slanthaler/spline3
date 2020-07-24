/**
   Unit tests for tridiagonal_matrix.hpp
**/

#pragma once
#include "gtest/gtest.h"
#include "tridiagonal_matrix.hpp"
#include "random.hpp"

#include <math.h>
#include <vector>

namespace{

	// only test for doubles for now
	class TridiagonalMatrixDouble : public ::testing::Test{ 
	protected:
		const int N=20;
		std::vector<double> l, d, u;
		std::vector<double> solution;
		TridiagonalMatrix<double> A;
	   
		TridiagonalMatrixDouble():
			l(N,1.0), d(N,4.0), u(N,-2.0)
		{
			A.resize(N);
			solution.resize(N);
			Random::fill_uniform(solution);
		};
	};

	TEST_F(TridiagonalMatrixDouble,SolveSymmetric){
		for(int i=0; i<N; ++i){
			A.upper[i] = l[i]; // symmetric! u==l
			A.diagonal[i] = d[i];
			A.lower[i] = l[i]; // symmetric!
		}
		auto b = A.dot(solution); // set up right-hand side
		A.solve_inplace(b);       // solve (should: b = solution)

		for(int i=0; i<N; ++i){
			EXPECT_NEAR(solution[i],b[i],1e-12);
		}
	};
	
	TEST_F(TridiagonalMatrixDouble,SolveAsymmetric){
		for(int i=0; i<N; ++i){
			A.upper[i] = u[i];
			A.diagonal[i] = d[i];
			A.lower[i] = l[i]; // asymmetric!
		}
		auto b = A.dot(solution); // set up right-hand side
		A.solve_inplace(b);       // solve (should: b = solution)

		for(int i=0; i<N; ++i){
			EXPECT_NEAR(solution[i],b[i],1e-12);
		}
	};

	TEST_F(TridiagonalMatrixDouble,SolveRandom){
		// fill with random (-1,1)-distributed
		Random::fill_uniform(u,-1,1);
		Random::fill_uniform(d,-1,1);
		Random::fill_uniform(l,-1,1);

		// make sure, we end up with solvable system
		for(int i=0; i<N; ++i){
			d[i] += 5.0;            // --> |d[i]| > 4.
			d[i] *= Random::sign(); // give a random sign
		};
		
		// fill matrix A
		for(int i=0; i<N; ++i){
			A.upper[i] = u[i];
			A.diagonal[i] = d[i];
			A.lower[i] = l[i];
		}
		auto b = A.dot(solution); // set up right-hand side
		A.solve_inplace(b);       // solve (should: b = solution)

		for(int i=0; i<N; ++i){
			EXPECT_NEAR(solution[i],b[i],1e-12);
		}
	};

} // namespace
