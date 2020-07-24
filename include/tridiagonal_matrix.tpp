/*
  tridiagonal_matrix.tpp
*/
#include <cassert>

template<typename T>
TridiagonalMatrix<T>::TridiagonalMatrix(int dim)
{
	resize(dim);
};

template<typename T>
void TridiagonalMatrix<T>::resize(int dim){
	dim_ = dim;
	upper.resize(dim);
	diagonal.resize(dim);
	lower.resize(dim);
};

template<typename T>
int TridiagonalMatrix<T>::dim() const{
	return dim_;
};

template<typename T>
std::vector<T> TridiagonalMatrix<T>::dot(const std::vector<T>& x){
	assert(x.size()==dim_);
	std::vector<T> result(dim_);
	int n = dim_-1;
	
	result[0] = diagonal[0]*x[0] + upper[0]*x[0+1];
	for(int i=1; i<n; ++i){
		result[i] = lower[i]*x[i-1] + diagonal[i]*x[i] + upper[i]*x[i+1];
	}
	result[n] = lower[n]*x[n-1] + diagonal[n]*x[n];

	return result;
}

template<typename T>
void TridiagonalMatrix<T>::solve_inplace(std::vector<T>& b){
	/* Thomas algorithm
	   [https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonalmatrix_algorithm]
	*/
	auto buffer_ = upper;
	int n = dim_ - 1;
	
	/*
	  c-->buffer_
	  b-->diagonal
	  d-->b
	  a-->lower
	*/
	buffer_[0] /= diagonal[0];
	b[0] /= diagonal[0];

	for(int i=1; i<n; i++){
		T factor = diagonal[i] - lower[i]*buffer_[i-1];
		buffer_[i] /= factor;
		b[i] = (b[i] - lower[i]*b[i-1]) / factor;
	}
	T factor = diagonal[n] - lower[n]*buffer_[n-1];
	b[n] = (b[n] - lower[n]*b[n-1]) / factor;

	for(int i=n-1; i>=0; i--){
		b[i] -= buffer_[i]*b[i+1];
	}
};


