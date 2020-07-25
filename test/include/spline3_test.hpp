/**
   Unit tests for spline3.hpp
**/

#pragma once
#include "gtest/gtest.h"
#include "spline3.hpp"
#include "regression.hpp"
#include "io.hpp"
#include "random.hpp"

#include <math.h>
#include <vector>

namespace{
	using array = std::vector<double>;

	double fun(double x){ return sin(5.*x); }
	double fun_p(double x){ return 5.*cos(5.*x); }
	
	class spline3Test : public ::testing::Test{ 
	protected:
		int N;
		array x,f;
		double fpL,fpR;
		Spline3::Spline3<double> spl;
		
		spline3Test(){
			N = 20;
			interpolate_on_grid(N,
								Spline3::BCtype::clamped,
								Spline3::BCtype::clamped);
		}

		void interpolate_on_grid(int Nin, Spline3::BCtype bcL, Spline3::BCtype bcR){
			// initialize x, f(x) = sin(x)
			N = Nin;
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
			Spline3::BC<double> bc(bcL,fpL,bcR,fpR);
			spl.SetSpline(x,f,bc);
			
		}

		double get_max_error(){
			double dx = 1e-5;
			double xx = x[0];
			double err = 0;
			
			while(xx < x[N-1]){
				err = std::max(err, std::abs(fun(xx)-spl(xx)));
				xx += dx;
			}
			return err;
		}

		double get_mean_error(){
			double dx = 1e-5;
			double xx = x[0];
			double err = 0;
			
			while(xx < x[N-1]){
				err += std::abs(fun(xx)-spl(xx))*dx;
				xx += dx;
			}
			return err/(x[N-1]-x[0]);
		}

		double convergence_rate(std::vector<int> Nvals,
								std::vector<double> error){
			std::vector<double> h;
			for(int i=0; i<Nvals.size(); ++i){
				h.push_back(1./(double) Nvals[i]);
				// take the log!
				h[i] = std::log(h[i]);
				error[i] = std::log(error[i]);
			}
			
			return regression::slope(h,error);
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


	TEST_F(spline3Test, NotAKnotBC){
		// interpolate with not-a-knot BC
		interpolate_on_grid(100,
							Spline3::BCtype::notaknot,
							Spline3::BCtype::notaknot);

		//
		ASSERT_NEAR(f[0],spl(x[0]),1e-10);
		ASSERT_NEAR(f[N-1],spl(x[N-1]),1e-10);
		ASSERT_NEAR(fpL,spl.p(x[0]),1e-3);
		ASSERT_NEAR(fpR,spl.p(x[N-1]),1e-3);
	};


	/*
	  Expected order of convergence is 4.
	 */
	TEST_F(spline3Test, ConvergenceTest){
		std::vector<int> Nvals;
		std::vector<double> max_error, mean_error;
		for(int N=10; N<5000; N*=2){
			Spline3::BCtype bcL = Spline3::BCtype::notaknot;
			Spline3::BCtype bcR = Spline3::BCtype::notaknot;
			interpolate_on_grid(N,bcL,bcR);
			//
			Nvals.push_back(N);
			max_error.push_back( get_max_error() );
			mean_error.push_back( get_mean_error() );
		}

		IO::write("max_error.dat", max_error);
		IO::write("mean_error.dat", mean_error);

		double pmax  = convergence_rate(Nvals,max_error);
		double pmean = convergence_rate(Nvals,mean_error);

	    EXPECT_GT(pmax,3.8);
	    EXPECT_GT(pmean,3.8);

		std::cout << "convergence rate: " << "pmax  = " << pmax << "\n";
		std::cout << "                  " << "pmean = " << pmean << "\n";
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

	// interpolation on non-equidistant grid
	TEST_F(spline3Test, ExactInterpolationNotEquidist) {
		Spline3::Spline3<double> spl;
		// initialize x, f(x) = sin(x)
		int N=0;
		double dx = 1e-1;
		double xnew = 0.;
		std::vector<double> x;
		std::vector<double> f;
		
		//
		x.push_back(0.);
		f.push_back(polynomial(0.0));
		++N;
		//
		while(xnew < 1.){
			xnew += dx*Random::uniform(0.1,1.0);
			if(xnew > 1.) xnew = 1.;
			x.push_back(xnew);
			f.push_back(polynomial(xnew));
			++N;
		}
		
		// save also boundary derivative values
		double fpL = polynomial_p(x[0]);
		double fpR = polynomial_p(x[N-1]);
		
		// interpolate
		Spline3::BC<double> bc(Spline3::BCtype::clamped,fpL,
							   Spline3::BCtype::clamped,fpR);
		spl.SetSpline(x,f,bc);
		spl.WriteOut("spline.dat");
		
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
