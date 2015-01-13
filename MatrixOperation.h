#pragma once
#include <flens/flens.cxx>

template<typename T, int N>
class MatrixOperation
{
public:
	template<typename Matrix>
	T Determinant(Matrix& M)
	{
		/*
		if (M.numRows() == 4)
		{
			return Determinant4x4(M);
		}
		else
		*/
		{
			typedef typename Matrix::IndexType IndexType;
			typedef flens::DenseVector< flens::Array<IndexType> > IndexVector;
			IndexVector piv(M.numRows());
			//Calculate LU factorization
			flens::lapack::trf(M, piv);
			T det = T(1.0);
			for (int i = M.firstCol(); i <= M.lastCol(); ++i)
				det *= M(i, i) * ((piv(i) != i) ? -T(1.0) : T(1.0));
			return det;
		}
	}
	
	template<typename Matrix, typename Pivot>
	T Determinant(Matrix& M, Pivot& piv)
	{
		//Calculate LU factorization
		flens::lapack::trf(M, piv);
		T det = T(1.0);
		for (int i = M.firstCol(); i <= M.lastCol(); ++i)
			det *= M(i, i) * ((piv(i) != i) ? -T(1.0) : T(1.0));
		return det;
	}
	
	//Computes inverse
	template<typename Matrix>
	void Inverse(Matrix& M)
	{
		typedef typename Matrix::IndexType IndexType;
		typedef flens::DenseVector< flens::Array<IndexType> > IndexVector;
		IndexVector piv(M.numRows());
		//Calculate LU factorization
		flens::lapack::trf(M, piv);
		//Use LU factorization to calculate inverse
		flens::lapack::tri(M, piv);
	}
	
	//Computes inverse from LU factorization as input
	template<typename Matrix, typename Pivot>
	void Inverse(Matrix& M, Pivot& piv)
	{
		//Use LU factorization to calculate inverse
		flens::lapack::tri(M, piv);
	}
public:
	template<typename Matrix>
	T Determinant4x4(Matrix& M)
	{
		return
			  M(1,1) * M(2,2) * M(3,3) * M(4,4) - M(1,1) * M(2,2) * M(3,4) * M(4,3) + M(1,1) * M(2,3) * M(3,4) * M(4,2)
			- M(1,1) * M(2,3) * M(3,2) * M(4,4) + M(1,1) * M(2,4) * M(3,2) * M(4,3) - M(1,1) * M(2,4) * M(3,3) * M(4,2)
			- M(1,2) * M(2,3) * M(3,4) * M(4,1) + M(1,2) * M(2,3) * M(3,1) * M(4,4) - M(1,2) * M(2,4) * M(3,1) * M(4,3)
			+ M(1,2) * M(2,4) * M(3,3) * M(4,1) - M(1,2) * M(2,1) * M(3,3) * M(4,4) + M(1,2) * M(2,1) * M(3,4) * M(4,3)
			+ M(1,3) * M(2,4) * M(3,1) * M(4,2) - M(1,3) * M(2,4) * M(3,2) * M(4,1) + M(1,3) * M(2,1) * M(3,2) * M(4,4)
			- M(1,3) * M(2,1) * M(3,4) * M(4,2) + M(1,3) * M(2,2) * M(3,4) * M(4,1) - M(1,3) * M(2,2) * M(3,1) * M(4,4)
			- M(1,4) * M(2,1) * M(3,2) * M(4,3) + M(1,4) * M(2,1) * M(3,3) * M(4,2) - M(1,4) * M(2,2) * M(3,3) * M(4,1)
			+ M(1,4) * M(2,2) * M(3,1) * M(4,3) - M(1,4) * M(2,3) * M(3,1) * M(4,2) + M(1,4) * M(2,3) * M(3,2) * M(4,1);
	}
};

template<typename T>
class MatrixOperation<T, 2>
{
public:
	template<typename Matrix, typename Pivot>
	T Determinant(Matrix& M, Pivot& piv)
	{
		det = M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1);
		validDeterminant = true;
		return det;
	}
	
	//Specialized version does require LU factorization as input
	//Uses determinant of last call to Determinant(M, piv)
	template<typename Matrix, typename Pivot>
	void Inverse(Matrix& M, Pivot& piv)
	{
		if (!validDeterminant)
			det = Determinant(M, piv);
		T a = M(1, 1);
		M(1, 1) = M(2, 2) / det;
		M(1, 2) = -M(1, 2) / det;
		M(2, 1) = -M(2, 1) / det;
		M(2, 2) = a / det;
	}
private:
	T det;
	bool validDeterminant = false;
};
