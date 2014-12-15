#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <list>
#include <utility>
#include <cmath>
#include <numeric>
#include <cstdint>
#include "LookUpTable.h"
//#define EIGEN_USE_MKL_ALL
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"
#include <flens/flens.cxx>
#include "measurements.h"
#include "random.h"
#include "parser.h"
#include "types.h"
#include "ConfigSpace.h"
#include "VertexHandler.h"

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
};

template<typename ConfigSpace_t>
class UpdateHandler
{
	public:
		typedef typename ConfigSpace_t::uint_t uint_t;
		typedef typename ConfigSpace_t::int_t int_t;
		typedef typename ConfigSpace_t::value_t value_t;
		typedef VertexHandler<ConfigSpace_t> VertexHandler_t;
		
		typedef flens::GeMatrix< flens::FullStorage<value_t, flens::ColMajor> > GeMatrix;
		typedef typename GeMatrix::IndexType IndexType;
		typedef flens::DenseVector< flens::Array<IndexType> > IndexVector;
		
		UpdateHandler(ConfigSpace_t& configSpace)
			: configSpace(configSpace), vertexHandler(VertexHandler_t(configSpace)), invG(GeMatrix(0, 0))
		{}
		
		void Init()
		{
			for (int_t n = 1; n <= 2; ++n)
				prefactor.push_back(std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), n));
		}
		
		//Performs an LU factorizatin on M
		value_t Determinant(GeMatrix& M, IndexVector& piv)
		{
			//Calculate LU factorization
			flens::lapack::trf(M, piv);
			value_t det = 1.0;
			for (int_t i = M.firstCol(); i <= M.lastCol(); ++i)
				det *= M(i, i) * ((piv(i) != i) ? -1.0 : 1.0);
			return det;
		}

		//Specialized version does not perform LU factorization
		value_t Determinant2x2(GeMatrix& M)
		{
			return M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1);
		}

		//Computes inverse from LU factorization as input
		void Inverse(GeMatrix& M, IndexVector& piv)
		{
			//Use LU factorization to calculate inverse
			flens::lapack::tri(M, piv);
		}

		//Specialized version does require LU factorization as input
		void Inverse2x2(GeMatrix& M, value_t det)
		{
			value_t a = M(1, 1);
			M(1, 1) = M(2, 2) / det;
			M(1, 2) = -M(1, 2) / det;
			M(2, 1) = -M(2, 1) / det;
			M(2, 2) = a / det;
		}

		template<int_t N>
		bool AddVertices()
		{
			uint_t k = vertexHandler.Vertices();
			GeMatrix u(2 * k, 2 * N);
			GeMatrix v(2 * N, 2 * k);
			GeMatrix a(2 * N, 2 * N);
			vertexHandler.template ComputeWoodburyMatrices<N>(u, v, a);

			GeMatrix invGu = invG * u;
			GeMatrix S = a - v * invGu;
			IndexVector pivS(2 * N);
			
			value_t detS = (N == 1 ? Determinant2x2(S) : Determinant(S, pivS));
			value_t preFactor = configSpace.AdditionFactorialRatio(k, N) * prefactor[N-1];
			value_t acceptRatio = preFactor * detS;
			if (acceptRatio < 0.0)
				std::cout << "AddVertices: AcceptRatio: " << acceptRatio << std::endl;
			if (configSpace.rng() < acceptRatio)
			{
				if(N == 1)
					Inverse2x2(S, detS);
				else
					Inverse(S, pivS);
				GeMatrix Q = -invGu * S;
				GeMatrix vinvG = v * invG;
				GeMatrix R = -S * vinvG;
				GeMatrix P = invG - invGu * R;
				
				invG.resize(2 * (k + N), 2 * (k + N));
				const flens::Underscore<IndexType> _;
				invG(_(1, 2 * k), _(1, 2 * k)) = P;
				invG(_(1, 2 * k), _(2 * k + 1, 2 * (k + N))) = Q;
				invG(_(2 * k + 1, 2 * (k + N)), _(1, 2 * k)) = R;
				invG(_(2 * k + 1, 2 * (k + N)), _(2 * k + 1, 2 * (k + N))) = S;
				vertexHandler.AddBufferedVertices();
				return true;
			}
			return false;
		}
		
		template<int_t N>
		bool RemoveVertices()
		{
			uint_t k = vertexHandler.Vertices();
			if (k < N)
				return false;
			
			GeMatrix perm(2 * k, 2 * k);
			vertexHandler.PermutationMatrix(perm);
			GeMatrix P = invG * perm;
			invG = flens::transpose(perm) * P;
			
			const flens::Underscore<IndexType> _;
			GeMatrix S = invG(_(2 * (k - N) + 1, 2 * k), _(2 * (k - N) + 1, 2 * k));
			IndexVector pivS(2 * N);
			value_t detS = (N == 1 ? Determinant2x2(S) : Determinant(S, pivS));
			value_t preFactor = configSpace.RemovalFactorialRatio(k, N) / prefactor[N-1];
			value_t acceptRatio = preFactor * detS;
			if (acceptRatio < 0.0)
				std::cout << "RemoveVertex: AcceptRatio" << acceptRatio << std::endl;
			if (configSpace.rng() < acceptRatio)
			{
				GeMatrix P = invG(_(1, 2 * (k - N)), _(1, 2 * (k - N)));
				GeMatrix Q = invG(_(1, 2 * (k - N)), _(2 * (k - N) + 1, 2 * k));
				GeMatrix R = invG(_(2 * (k - N) + 1, 2 * k), _(1, 2 * (k - N)));
				if(N == 1)
					Inverse2x2(S);
				else
					Inverse(S, pivS);
				
				invG.resize(2 * (k - N), 2 * (k - N));
				GeMatrix SR = S * R;
				invG = P - Q * SR;
				
				vertexHandler.RemoveBufferedVertices();
				return true;
			}
			else
			{
				GeMatrix O = perm * invG;
				invG = O * flens::transpose(perm);
				return false;
			}
		}
		
		void SymmetrizeInvG()
		{
			SymmetrizeMatrix(invG);
		}
		
		void SymmetrizeMatrix(GeMatrix& M)
		{	
			for (uint_t i = M.firstRow(); i <= M.lastRow(); ++i)
			{
				for (uint_t j = M.firstCol(); j <= i; ++j)
				{
					value_t mean = (std::abs(M(i, j)) + std::abs(M(j, i))) / 2.0;
					M(i, j) = sgn(M(i, j)) * mean;
					M(j, i) = sgn(M(j, i)) * mean;
				}
				M(i, i) = 0.;
			}
		}
		
		VertexHandler_t& GetVertexHandler()
		{
			return vertexHandler;
		}
		
	private:
		ConfigSpace_t& configSpace;
		VertexHandler_t vertexHandler;
		GeMatrix invG;
		std::vector<value_t> prefactor;
};