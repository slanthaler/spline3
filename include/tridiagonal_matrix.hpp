#pragma once
/*
  tridiagonal_matrix.hpp
*/
#include <vector>
#include <iostream>
#include <iomanip>

template<typename T>
class TridiagonalMatrix
{
public:
    TridiagonalMatrix() {};     // constructor
    TridiagonalMatrix(int dim); // constructor
    ~TridiagonalMatrix() {};    // destructor
    void resize(int dim);       // init with dim,n_u,n_l
    int dim() const;            // matrix dimension
    // 
    //std::vector<T> solve(const std::vector<T>& b);
	void solve_inplace(std::vector<T>& b);
	std::vector<T> dot(const std::vector<T>& x);

	//
	void print_matrix() const{
		for(int i=0; i<dim_; ++i){
			for(int j=0; j<dim_; ++j){
				std::cout << std::setw(6) << std::setprecision(2) << std::fixed;
				if(j==i-1)
					std::cout << lower[i];
				else if(j==i)
					std::cout << diagonal[i];
				else if(j==i+1)
					std::cout << upper[i];
				else
					std::cout <<  0.;
				//
				std::cout << ", ";
			}
			std::cout << "\n";
		}

	}

	//
	void print_matrix(const std::vector<T>& b) const{
		for(int i=0; i<dim_; ++i){
			for(int j=0; j<dim_; ++j){
				std::cout << std::setw(6) << std::setprecision(2) << std::fixed;
				if(j==i-1)
					std::cout << lower[i];
				else if(j==i)
					std::cout << diagonal[i];
				else if(j==i+1)
					std::cout << upper[i];
				else
					std::cout <<  0.;
				//
				std::cout << ", ";
			}

			//
			std::cout << "   |   " << b[i];
			
			std::cout << "\n";
		}
		
	}

private:
	int dim_;
public: // TODO
	std::vector<T> upper;
	std::vector<T> diagonal;
	std::vector<T> lower;
};

#include "tridiagonal_matrix.tpp"
