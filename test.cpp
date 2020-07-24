#include <iostream>
#include "tridiagonal_matrix.hpp"

int main(){
	const int dim = 4;
	TridiagonalMatrix<double> A(dim);

	/*
	  Corresponding to the matrix
	  
	  A = {  4, -2,  0,  0
	        -1,  4, -2,  0
		 	 0, -1,  4, -2
			 0,  0,  -1, 4 }

	  b = {  2,
	         0, 
		     3, 
		    23  }

	  x = {  2,
	         3,
			 5,
			 7  }
	 */
    std::vector<double> l { 0, -1, -1, -1 };
	std::vector<double> d { 4,  4,  4,  4 };
	std::vector<double> u {-2, -2, -2,  0 };
	std::vector<double> b { 2,  0,  3, 23 };
	//
	std::vector<double> result { 2,  3,  5, 7  };
	
	for(int i=0; i<dim; ++i){
		A.lower(i)    = l[i];
		A.diagonal(i) = d[i];
		A.upper(i)    = u[i];
	}

	auto b_check = A.dot(result);
	for(int i=0; i<dim; ++i){
		std::cout << "A.x: " << b_check[i]
				  << "  ---: " << b[i]
				  << "\n";
	}

	// solve the system
	A.solve_inplace(b);

	for(int i=0; i<dim; ++i){
		std::cout << "computed: " << b[i]
				  << "  -- actual: " << result[i] << "\n";
	}
	
	return 0;
}
