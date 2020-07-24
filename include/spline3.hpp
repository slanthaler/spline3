/**
   spline3.hpp

   Implements cubic Spline interpolation class.
   The splines are locally represented on an interval [x_i,x_{i+1}) 
   in the form:
   $$ y(x) = a_i + b_i (x-x_i) + c_i (x-x_i)^2 + d_i (x-x_i)^3.$$
   The coefficients are determined from 

   1.) a_i=y_i (given by input),
   2.) y'(x_i-)  = y'(x_i+)  for all i=0,...,N,
   3.) y''(x_i-) = y''(x_i+) for all i=0,...,N.

**/

#pragma once
#include<vector>
#include<string>
#include<fstream>

#include "tridiagonal_matrix.hpp"

using Grid = std::vector<double>;
template<typename T>
using Values = std::vector<T>;

namespace Spline3{
	/* 
	   struct to define commonly used boundary condition types
	 */
    enum BCtype{notaknot,clamped};
	template<typename T> struct BC{
		BCtype left_type, right_type;
		T left_val, right_val;

		BC(BCtype left_t,  T left_v,
		   BCtype right_t, T right_v):
			left_type(left_t),
			left_val(left_v),
			right_type(right_t),
			right_val(right_v)
		{};
		
	};

	/*
	  struct to store data for tridiagonal system
	 */
	template<typename T>
	struct TridiagonalSystem{
		TridiagonalMatrix<T> A;
		std::vector<T> b, h, s;
	};
	
	/*
	  Cubic spline interpolation class
	 */
	template<typename T>
	class Spline3{

	public:
		/* 
		   Constructor 
		*/
		Spline3(){};
		Spline3(Grid grid, Values<T> vals);

		/*
		  After-the-fact-constructor
		*/
		void SetSpline(Grid grid, Values<T> vals, BC<T> bc);
		
		/* 
		   Overload operator() for evaluation at x
		*/
		T operator()(double x) const;
		Values<T> operator()(Grid x) const;

		/*
		  Find (left) index corresponding to interval
		 */
		int index(double x) const{
			int iL = 0;
			int iR = grid_.size()-1;
			if(x < grid_[iL])  return iL;
			if(x >= grid_[iR]) return iR-1;
			while(iR-iL > 1){
				int imid = (iR+iL)/2;
				if(x < grid_[imid])
					iR = imid;
				else
					iL = imid;
			}
			return iL;
		};
		
		/* 
		   Evaluate derivative [p] and second derivative [pp] at x
		*/
		// single-point evaluation
		T p(double x) const;    
		T pp(double x) const;
		// evaluate on "grid" of points
		Values<T> p(Grid x) const; 
		Values<T> pp(Grid x) const; 

		/*
		  Evaluate integral [I] at x (initialized to 0 at left boundary)
		*/
		T I(double x) const;
		Values<T> I(Grid x) const;
		
		/*
		  Solution of the spline interpolation system
		 */
		void SetUpSystem(TridiagonalSystem<T>& system);
		void ApplyBC(TridiagonalSystem<T>& system, BC<T> bc);
		void SolveSystem(TridiagonalSystem<T>& system);
		void GetCoefficients(const TridiagonalSystem<T>& system);

		/*
		  Write to file
		 */
		void WriteOut(std::string filename){
			if(grid_.size()<1) return;

			std::ofstream file;
			file.open(filename);
			for(int i=0; i<grid_.size()-1;++i){
				int Nfine = 20;
				double dx = (grid_[i+1]-grid_[i])/Nfine;
				for(int j=0; j<Nfine; ++j){
					double x = grid_[i] + j*dx;
					file << x
						 << ", " << (*this)(x)
						 << ", " << (*this).p(x)
						 << ", " << (*this).pp(x)
						 << "\n";
				}
			}
			file.close();
		}

		/*
		  Write coefficients to file
		 */
		void WriteOutCoefficients(std::string filename){
			if(grid_.size()<1) return;

			std::ofstream file;
			file.open(filename);
			for(int i=0; i<grid_.size()-1;++i){			
				file << i << ", "
					 << coeffs_[i+0] << ", "
					 << coeffs_[i+1] << ", "
					 << coeffs_[i+2] << ", "
					 << coeffs_[i+3] << ", " 
					 << "\n";
			}
			file.close();
		}

		
	private:
		std::vector<T> grid_;
		std::vector<T> vals_;
		std::vector<T> coeffs_;
		std::vector<T> int_;
	};


} // namespace


#include "spline3.tpp"
