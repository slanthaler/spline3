/**
   Unit tests for spline3.hpp
**/

#pragma once
#include "gtest/gtest.h"
#include "spline3.hpp"

#include <math.h>
#include <vector>

namespace{
	using array = std::vector<double>;

	double fun(double x){ return sin(5.*x); }
	double fun_p(double x){ return 5.*cos(5.*x); }
	
	class spline3Test : public ::testing::Test{ 
	protected:
		const int N=20;
		array x,f;
		double fpL,fpR;
		Spline3::Spline3<double> spl;
		
		spline3Test(){
			// initialize x, f(x) = sin(x)
			double dx = 1./double(N-1);
			x.resize(N);
			f.resize(N);
			for(int i=0; i<N; ++i){
				x[i] = (i-1)*dx;
				f[i] = fun(x[i]);
			}
			
			// save also boundary derivative values
		    fpL = fun_p(x[0]);
			fpR = fun_p(x[N-1]);

			// interpolate
			Spline3::BC<double> bc(Spline3::BCtype::clamped,fpL,
								   Spline3::BCtype::clamped,fpR);
			spl.SetSpline(x,f,bc);

			// FIXME write to file
			spl.WriteOut("spline.dat");
			spl.WriteOutCoefficients("coeffs.dat");
		}
	};


	TEST_F(spline3Test, indexMatching) {
		for(int i=0; i<N-1; ++i){
			EXPECT_EQ(i,spl.index(x[i]));
		}
	};
	
	TEST_F(spline3Test, indexToLeft) {
		for(int i=0; i<N-1; ++i){
			double xmid = (x[i] + x[i+1])/2.;
			EXPECT_EQ(i,spl.index(xmid));
		}
	};

	TEST_F(spline3Test, indexOutOfBounds) {
		EXPECT_EQ(0,spl.index(x[0]-1.));
		EXPECT_EQ(N-2,spl.index(x[N-1]+1.)); // last interval is N-2!
		EXPECT_EQ(N-2,spl.index(x[N-1]));
	};
	
	TEST_F(spline3Test, InterpolatedValues) {
		for(int i=0; i<N; ++i){
		    ASSERT_NEAR(spl(x[i]),f[i],1e-10);
		}
	};

	TEST_F(spline3Test, Continuity) {
		double tiny=1e-13;
		for(int i=0; i<N; ++i){
			ASSERT_NEAR(spl(x[i]+tiny),
						spl(x[i]-tiny),
						1e-10);
		}
	};

	TEST_F(spline3Test, C1Continuity) {
		double tiny=1e-13;
		for(int i=0; i<N; ++i){
			ASSERT_NEAR(spl.p(x[i]+tiny),
						spl.p(x[i]-tiny),
						1e-10);
		}
	};

	TEST_F(spline3Test, C2Continuity) {
		double tiny=1e-13;
		for(int i=0; i<N; ++i){
			ASSERT_NEAR(spl.pp(x[i]+tiny),
						spl.pp(x[i]-tiny),
						1e-10);
		}
	};

	TEST_F(spline3Test, BoundaryConditions) {
		ASSERT_NEAR(fpL,spl.p(x[0]),1e-12);
	    ASSERT_NEAR(fpR,spl.p(x[N-1]),1e-12);
	};



	
	inline double polynomial(double x)
	{
		return 1. + x + 2.*x*x;
	};

	inline double polynomial_p(double x)
	{
		return 1. + 4.*x;
	};
	
	TEST_F(spline3Test, ExactInterpolation) {
		Spline3::Spline3<double> spl;
		// initialize x, f(x) = sin(x)
		int N=20;
		double dx = 1./double(N-1);
		array x(N);
		array f(N);
		for(int i=0; i<N; ++i){
			x[i] = (i-1)*dx;
			f[i] = polynomial(x[i]);
		}
		
		// save also boundary derivative values
		double fpL = polynomial_p(x[0]);
		double fpR = polynomial_p(x[N-1]);
		
		// interpolate
		Spline3::BC<double> bc(Spline3::BCtype::clamped,fpL,
							   Spline3::BCtype::clamped,fpR);
		spl.SetSpline(x,f,bc);

		//
		dx /= 100.;
		double x_f = x[0];
		while(x_f < x[N-1]){
			ASSERT_NEAR(polynomial(x_f),
						spl(x_f),
						1e-12);
			x_f += dx;
		};
	};

} // namespace
