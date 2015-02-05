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
		template<int_t N, int_t M> using matrix_t = Eigen::Matrix<value_t, N, M, Eigen::RowMajor>;
		template<int_t N> using inv_solver_t = Eigen::FullPivLU< matrix_t<N, N> >;
		
		struct Matrices
		{
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> invG;
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormU;
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormV;
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormA;
			value_t detWormS;
		};
		
		UpdateHandler(ConfigSpace_t& configSpace)
			: configSpace(configSpace), vertexHandler(VertexHandler_t(configSpace))
		{}
		
		void Init()
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
			std::cout << M.format(CleanFmt) << std::endl;
			//std::cout << M.format(OctaveFmt) << std::endl;
			//std::cout << M.format(HeavyFmt) << std::endl;
		}

		template<int_t N>
		bool AddVertices(value_t preFactor, bool isWorm)
		{
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;
			matrix_t<Eigen::Dynamic, n> u(k, n);
			matrix_t<n, Eigen::Dynamic> v(n, k);
			matrix_t<n, n> a(n, n);
			vertexHandler.WoodburyAddVertices(u, v, a);

			matrix_t<Eigen::Dynamic, n> invGu = invG * u;
			matrix_t<n, n> invS = a;
			invS.noalias() -= v * invGu;
	
			value_t acceptRatio = preFactor * invS.determinant();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "AddVertices(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
			}
			if (configSpace.rng() < acceptRatio)
			{
				matrix_t<n, n> S = invS.inverse();
				matrix_t<n, Eigen::Dynamic> R = -S * v * invG;
				
				invG.conservativeResize(k + n, k + n);
				invG.topLeftCorner(k, k).noalias() -= invGu * R;
				invG.topRightCorner(k, n).noalias() = -invGu * S;
				invG.bottomLeftCorner(n, k) = R;
				invG.template bottomRightCorner<n, n>() = S;
				
				vertexHandler.AddBufferedVertices(isWorm);
				return true;
			}
			return false;
		}
		
		template<int_t N>
		bool RemoveVertices(double preFactor, bool isWorm)
		{
			if (isWorm && vertexHandler.Worms() < N)
				return false;
			if ((!isWorm) && vertexHandler.Vertices() < N)
				return false;
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;

			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(k);
			vertexHandler.PermutationMatrix(perm.indices(), isWorm);
			invG = perm.transpose() * invG * perm;
			
			matrix_t<n, n> S = invG.template bottomRightCorner<n, n>();
			value_t acceptRatio = preFactor * S.determinant();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "RemoveVertices(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
				vertexHandler.PrintVertices();
				vertexHandler.PrintWormVertices();
				vertexHandler.PrintIndexBuffer();
				for (uint_t i = 0; i < k; ++i)
					std::cout << perm.indices()[i] << " ";
				std::cout << std::endl;
				std::cout << preFactor << std::endl;
				std::cout << S.determinant() << std::endl;
				std::cout << "---" << std::endl;
			}
			if (configSpace.rng() < acceptRatio)
			{
				matrix_t<n, n> invS = S.inverse();
				invG.topLeftCorner(k - n, k - n).noalias() -= invG.topRightCorner(k - n, n) * invS * invG.bottomLeftCorner(n, k - n);
				invG.conservativeResize(k - n, k - n);

				vertexHandler.RemoveBufferedVertices(isWorm);
				return true;
			}
			else
			{
				invG = perm * invG * perm.transpose();
				return false;
			}
		}
		
		template<int_t W>
		bool ShiftWorm()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t l = 2 * W;

			matrix_t<Eigen::Dynamic, l> wormU(k, l), shiftedWormU(k, l);
			matrix_t<l, Eigen::Dynamic> wormV(l, k), shiftedWormV(l, k);
			matrix_t<l, l> wormA(l, l), shiftedWormA(l, l);

			vertexHandler.ShiftWormToBuffer();
			vertexHandler.WoodburyWorm(wormU, wormV, wormA);
			vertexHandler.WoodburyShiftWorm(shiftedWormU, shiftedWormV, shiftedWormA);

			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(k + l);
			vertexHandler.PermutationMatrix(perm.indices(), true);
			invG = perm.transpose() * invG * perm;
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> M = invG.topLeftCorner(k, k);
			M.noalias() -= invG.topRightCorner(k, l) * invG.template bottomRightCorner<l, l>().inverse() * invG.bottomLeftCorner(l, k);

			matrix_t<l, l> invS = wormA;
			invS.noalias() -= wormV * M * wormU;
			value_t detInvS = invS.determinant();

			matrix_t<l, l> shiftedInvS = shiftedWormA;
			shiftedInvS.noalias() -= shiftedWormV * M * shiftedWormU;
			value_t detShiftedInvS = shiftedInvS.determinant();

			value_t acceptRatio = detShiftedInvS / detInvS * vertexHandler.WormShiftParity();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "WormShift: AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "Worms:" << std::endl;
				vertexHandler.PrintWormVertices();
			}
			if (print && configSpace.rng() < acceptRatio)
			{
				matrix_t<l, l> S = shiftedInvS.inverse();
				matrix_t<l, Eigen::Dynamic> R = -S * shiftedWormV * M;
				matrix_t<Eigen::Dynamic, l> Mu = M * shiftedWormU;
				invG.topLeftCorner(k, k) = M;
				invG.topLeftCorner(k, k).noalias() -= Mu * R;
				invG.topRightCorner(k, l).noalias() = -Mu * S;
				invG.bottomLeftCorner(l, k) = R;
				invG.template bottomRightCorner<l, l>() = S;
				invG = perm * invG * perm.transpose();

				vertexHandler.ApplyWormShift();
				return true;
			}
			else
			{
				invG = perm * invG * perm.transpose();
				return false;
			}
		}
		
		template<typename Matrix>
		value_t MatrixCondition(Matrix& M)
		{
			if (M.rows() == 0)
				return 0.0;
			Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
			std::cout << svd.singularValues()(M.rows()-1) << std::endl;
			return svd.singularValues()(0) / svd.singularValues()(M.rows()-1);
		}

		void SymmetrizeInvG()
		{
			SymmetrizeMatrix(invG);
		}

		value_t StabilizeInvG()
		{
			if (invG.rows() == 0)
				return 0.0;
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> G(invG.rows(), invG.cols());
			vertexHandler.PropagatorMatrix(G);
			inv_solver_t<Eigen::Dynamic> solver(G);
			invG = solver.inverse();
			return 0.0;
		}
		
		value_t StabilizeInvG(value_t& avgError, value_t& relError)
		{
			if (invG.rows() == 0)
				return 0.0;
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> G(invG.rows(), invG.cols());
			vertexHandler.PropagatorMatrix(G);
			inv_solver_t<Eigen::Dynamic> solver(G);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> stabInvG = solver.inverse();
			/*
			Eigen::JacobiSVD< matrix_t<Eigen::Dynamic, Eigen::Dynamic> > svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> sv = svd.singularValues().asDiagonal();
			for (uint_t i = 0; i < sv.rows(); ++i)
			{
				if (sv(i, i) < 0.00001)
					sv(i, i) = 0.0;
				else
					sv(i, i) = 1.0 / sv(i, i);
			}
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> stabInvG = svd.matrixV() * sv * svd.matrixU().transpose();
			*/
			avgError = 0.0;
			relError = 0.0;
			value_t N = stabInvG.rows() * stabInvG.rows();
			for (uint_t i = 0; i < stabInvG.rows(); ++i)
			{
				for (uint_t j = 0; j < stabInvG.cols(); ++j)
				{
					value_t err = std::abs(invG(i, j) - stabInvG(i, j));
					avgError += err / N;
					relError += std::abs(stabInvG(i, j)) / N;
				}
				relError = avgError / relError;
			}

			Eigen::JacobiSVD< matrix_t<Eigen::Dynamic, Eigen::Dynamic> > svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> sv = svd.singularValues();
			invG = stabInvG;
			return sv(0) / sv(G.rows()-1);
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
		
		void Serialize(odump& d)
		{
			vertexHandler.Serialize(d);
		}
		
		void Serialize(idump& d)
		{
			vertexHandler.Serialize(d);
			invG.resize(2 * vertexHandler.Vertices(), 2 * vertexHandler.Vertices());
			StabilizeInvG();
		}
		
	private:
		ConfigSpace_t& configSpace;
		VertexHandler_t vertexHandler;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> invG;
		uint_t maxWorms = 2;
		bool print = true;
};