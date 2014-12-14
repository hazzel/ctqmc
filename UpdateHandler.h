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
		template<int_t N, int_t M> using matrix_t = Eigen::Matrix<value_t, N, M>;
		
		UpdateHandler(ConfigSpace_t& configSpace)
			: configSpace(configSpace), vertexHandler(VertexHandler_t(configSpace))
		{}
		
		template<typename Matrix_t>
		void PrintMatrix(const Matrix_t& M)
		{
			using namespace Eigen;
			IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
			IOFormat CleanFmt(FullPrecision, 0, ", ", "\n", "[", "]");
			IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
			IOFormat HeavyFmt(FullPrecision, 0, ", ", ",\n", "{", "}", "{", "}");
			//std::cout << M.format(CommaInitFmt) << std::endl;
			//std::cout << M.format(CleanFmt) << std::endl;
			//std::cout << M.format(OctaveFmt) << std::endl;
			std::cout << M.format(HeavyFmt) << std::endl;
		}
		
		template<int_t N>
		bool AddVertices()
		{
			uint_t k = vertexHandler.Vertices();
			matrix_t<Eigen::Dynamic, 2 * N> u(2 * k, 2 * N);
			matrix_t<2 * N, Eigen::Dynamic> v(2 * N, 2 * k);
			matrix_t<2 * N, 2 * N> a(2 * N, 2 * N);
			vertexHandler.template ComputeWoodburyMatrices<N>(u, v, a);

			matrix_t<2 * N, 2 * N> invS;
			matrix_t<Eigen::Dynamic, 2 * N> invGu = invG * u;
			invS.noalias() = a - v * invGu;
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), N) * configSpace.AdditionFactorialRatio(k, N);
			value_t acceptRatio = preFactor * invS.determinant();
			if (acceptRatio < 0.0)
				std::cout << "AddVertices: AcceptRatio: " << acceptRatio << std::endl;
			if (configSpace.rng() < acceptRatio)
			{
				matrix_t<2 * N, 2 * N> S = invS.inverse();
				matrix_t<Eigen::Dynamic, 2 * N> Q;
				Q.noalias() = -invGu * S;
				matrix_t<2 * N, Eigen::Dynamic> R;
				R.noalias() = -S * v * invG;
				matrix_t<Eigen::Dynamic, Eigen::Dynamic> P;
				P.noalias() = invG - invGu * R;
				
				invG.resize(2 * (k + N), 2 * (k + N));
				invG.topLeftCorner(2 * k, 2 * k) = P;
				invG.topRightCorner(2 * k, 2 * N) = Q;
				invG.bottomLeftCorner(2 * N, 2 * k) = R;
				invG.template bottomRightCorner<2 * N, 2 * N>() = S;
				
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
			
			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm = vertexHandler.PermutationMatrix();
			invG = perm.transpose() * invG * perm;
			
			matrix_t<2 * N, 2 * N> S = invG.template bottomRightCorner<2 * N, 2 * N>();
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), -N) * configSpace.RemovalFactorialRatio(k, N);
			value_t acceptRatio = preFactor * S.determinant();
			if (acceptRatio < 0.0)
				std::cout << "RemoveVertex: AcceptRatio" << acceptRatio << std::endl;
			if (configSpace.rng() < acceptRatio)
			{
				matrix_t<Eigen::Dynamic, Eigen::Dynamic> P = invG.topLeftCorner(2 * (k - N), 2 * (k - N));
				matrix_t<Eigen::Dynamic, 2 * N> Q = invG.topRightCorner(2 * (k - N), 2 * N);
				matrix_t<2 * N, Eigen::Dynamic> R = invG.bottomLeftCorner(2 * N, 2 * (k - N));
				invG.resize(2 * (k - N), 2 * (k - N));
				invG.noalias() = P - Q * S.inverse() * R;
				
				vertexHandler.RemoveBufferedVertices();
				return true;
			}
			else
			{
				invG = perm * invG * perm.transpose();
				return false;
			}
		}
		
		void SymmetrizeInvG()
		{
			SymmetrizeMatrix(invG);
		}
		
		template<typename Matrix>
		void SymmetrizeMatrix(Matrix& M)
		{	
			for (uint_t i = 0; i < M.rows(); ++i)
			{
				for (uint_t j = 0; j < i; ++j)
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
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> invG;
		Eigen::FullPivLU< matrix_t<Eigen::Dynamic, Eigen::Dynamic> > invSolver;
};