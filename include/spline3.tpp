/* 
   spline3.tpp

   Implementation of cubic spline interpolation routines defined in spline3.hpp.
   Note: The Spline3 class is templated (allowing real/complex). Therefore, we
         cannot simply split these routines in .hpp/.cpp files, as usual. This 
         file must be included in spline3.hpp.
   Alternative: We could also have used a .cpp file, and then instantiated the 
         required instantiations in the .cpp; e.g. have
		     
		     template class Spline3<real_t>; 
		     template class Spline3<complex_t>;
		 
		 somehwere. But the current solution seems better to me!
 */
#include <iostream>

namespace Spline3{
	template<typename T>
	Spline3<T>::Spline3(Grid grid, Values<T> vals){
		SetSpline(grid,vals);
	};

	template<typename T>
	void Spline3<T>::SetSpline(Grid grid, Values<T> vals, BC<T> bc){
		// allocation
		grid_ = grid;
		vals_ = vals;
		coeffs_.resize(4*grid.size());

		//
		TridiagonalSystem<T> system;
	    SetUpSystem(system);
		ApplyBC(system,bc);
		SolveSystem(system);
	    GetCoefficients(system);
	};

	template<typename T>
	void Spline3<T>::SetUpSystem(TridiagonalSystem<T>& system){
		// allocation
		system.b.resize(grid_.size());
		system.h.resize(grid_.size());
		system.s.resize(grid_.size());
		system.A.resize(grid_.size());
		
		//
		for(int i=0; i<grid_.size()-1; ++i){
			system.h[i] = grid_[i+1] - grid_[i];
			system.s[i] = (vals_[i+1] - vals_[i])/system.h[i];
		}

		//
		system.A.lower[0] = 0.;
		system.A.diagonal[0] = 0.;
		system.A.upper[0] = 0.;
		//
		int N = grid_.size()-1;
		system.A.lower[N] = 0.;
		system.A.diagonal[N] = 0.;
		system.A.upper[N] = 0.;

		//
		for(int i=1; i<grid_.size()-1; ++i){
			// i-th equation
			system.A.lower[i]    = 1./system.h[i-1];
			system.A.upper[i]    = 1./system.h[i];
			system.A.diagonal[i] = 2.0*(1./system.h[i-1] + 1./system.h[i]);
			system.b[i] = 3.*(system.s[i-1]*system.A.lower[i]
							  + system.s[i]*system.A.upper[i]);
		}
	};

	template<typename T>
	void Spline3<T>::ApplyBC(TridiagonalSystem<T>& system,
							  BC<T> bc                    )
	{
		int N = system.A.dim()-1; // last available index
		switch(bc.left_type)
		{
		case BCtype::clamped:
			//system.A.lower[0] = 0.;
			system.A.diagonal[0] = 1.;
			system.A.upper[0] = 0.;
			system.b[0]   = bc.left_val;
		case BCtype::notaknot:
			//TODO;
			;
	    default:
			//TODO [error]
			;
		}
		switch(bc.right_type)
		{
		case BCtype::clamped:
			system.A.lower[N]    = 0.;
			system.A.diagonal[N] = 1.;
			system.b[N]          = bc.right_val;
		case BCtype::notaknot:
			//TODO
			;
		default:
			//TODO [error]
			;
		}
	};

	template<typename T>
	void Spline3<T>::SolveSystem(TridiagonalSystem<T>& system){
		system.A.solve_inplace(system.b);
	};

	template<typename T>
	void Spline3<T>::GetCoefficients(const TridiagonalSystem<T>& system){
		/*
		  Local representation of the form a_i+b_i*h+c_i*h^2+d_i*h^3.
		*/
		for(int i=0; i<system.b.size()-1; ++i){
			auto bi = system.b[i];
			auto bi1 = system.b[i+1];
			auto si = system.s[i];
			auto hi = system.h[i];
			//
			coeffs_[4*i+0] = vals_[i];                     // a_i
			coeffs_[4*i+1] = bi;                           // b_i
			coeffs_[4*i+2] = (-2*bi - bi1 + 3*si)/hi;      // c_i
			coeffs_[4*i+3] = (bi + bi1 - 2.*si) / (hi*hi); // d_i
		};
	};
	
	template<typename T>
	T Spline3<T>::operator()(double x) const{
		int i = index(x);
		double dx = x - grid_[i];
		int a = 4*i + 0;
		int b = 4*i + 1;
		int c = 4*i + 2;
		int d = 4*i + 3;
		return coeffs_[a] + dx*(coeffs_[b] + dx*(coeffs_[c] + dx*coeffs_[d]));
	};

	template<typename T>
	Values<T> Spline3<T>::operator()(Grid x) const{
		Values<T> result(x.size());
		for(int i=0; i<x.size(); ++i){
			result[i] = (*this)(x[i]);
		}
		// TODO: optimize me
	};

	template<typename T>
	T Spline3<T>::p(double x) const{
		int i = index(x);
		double dx = x - grid_[i];
		int b = 4*i + 1;
		int c = 4*i + 2;
		int d = 4*i + 3;
		return coeffs_[b] + dx*(2.*coeffs_[c] + 3*dx*coeffs_[d]);
	};    

	template<typename T>
	Values<T> Spline3<T>::p(Grid x) const{
		Values<T> result(x.size());
		for(int i=0; i<x.size(); ++i){
			result[i] = (*this).p(x[i]);
		}
		// TODO: optimize me
	}; 

	template<typename T>
	T Spline3<T>::pp(double x) const{
		int i = index(x);
		double dx = x - grid_[i];
		int c = 4*i + 2;
		int d = 4*i + 3;
		return 2.*coeffs_[c] + 6.*dx*coeffs_[d];
	};
	
	template<typename T>
	Values<T> Spline3<T>::pp(Grid x) const{
		Values<T> result(x.size());
		for(int i=0; i<x.size(); ++i){
			result[i] = (*this).pp(x[i]);
		}
		// TODO: optimize me
	};

	template<typename T>
	T Spline3<T>::I(double x) const{
		// TODO
	};
	
	template<typename T>
	Values<T> Spline3<T>::I(Grid x) const{
		// TODO
	};


} // namespace
