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
		bool AddVertices()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t n = 2 * N;
			matrix_t<Eigen::Dynamic, n> u(k, n);
			matrix_t<n, Eigen::Dynamic> v(n, k);
			matrix_t<n, n> a(n, n);
			vertexHandler.WoodburyAddVertices(u, v, a);

			matrix_t<Eigen::Dynamic, n> invGu = invG * u;
			matrix_t<n, n> invS = a;
			invS.noalias() -= v * invGu;
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), N) * configSpace.AdditionFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor * invS.determinant();
			if (acceptRatio < 0.0)
			{
				std::cout << "AddVertices(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "VertexBuffer:" << std::endl;
				vertexHandler.PrintVertexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				invSolver.compute(invS);
				matrix_t<n, n> S = invSolver.inverse();
				matrix_t<n, Eigen::Dynamic> R = -S * v * invG;
				
				invG.conservativeResize(k + n, k + n);
				invG.topLeftCorner(k, k).noalias() -= invGu * R;
				invG.topRightCorner(k, n).noalias() = -invGu * S;
				invG.bottomLeftCorner(n, k) = R;
				invG.template bottomRightCorner<n, n>() = S;
				
				vertexHandler.AddBufferedVertices();
				return true;
			}
			return false;
		}
		
		template<int_t N>
		bool AddVerticesWithWorms()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			uint_t l = 2 * vertexHandler.Worms();
			const uint_t n = 2 * N;
			
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> u(k, l + n);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> v(l + n, k);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> a(l + n, l + n);
			vertexHandler.WoodburyAddVertices(u, v, a);
			u.rightCols(l) = wormU;
			v.bottomRows(l) = wormV;
			a.bottomRightCorner(l, l) = wormA;

			matrix_t<Eigen::Dynamic, Eigen::Dynamic> invS = a;
			invS.noalias() -= v * invG * u;
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), N) * configSpace.AdditionFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor * invS.determinant() * detWormS;
			if (acceptRatio < 0.0)
			{
				std::cout << "AddVerticesWithWorm(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "VertexBuffer:" << std::endl;
				vertexHandler.PrintVertexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				wormU.conservativeResize(k + n, l);
				wormU.bottomRows(n) = a.topRightCorner(n, l);
				wormV.conservativeResize(l, k + n);
				wormV.rightCols(n) = a.bottomLeftCorner(l, n);

				u.conservativeResize(k, n);
				v.conservativeResize(n, k);
				a.conservativeResize(n, n);
				invS.resize(n, n);
				invS = a;
				invS.noalias() -= v * invG * u;
				
				invSolver.compute(invS);
				matrix_t<n, n> S = invSolver.inverse();
				matrix_t<Eigen::Dynamic, n> invGu = invG * u;
				matrix_t<n, Eigen::Dynamic> R = -S * v * invG;
				
				invG.conservativeResize(k + n, k + n);
				invG.topLeftCorner(k, k).noalias() -= invGu * R;
				invG.topRightCorner(k, n).noalias() = -invGu * S;
				invG.bottomLeftCorner(n, k) = R;
				invG.template bottomRightCorner<n, n>() = S;

				detWormS = 1.0 / (wormA - wormV * invG * wormU).determinant();
				vertexHandler.AddBufferedVertices();
				return true;
			}
			return false;
		}
		
		template<int_t N>
		bool RemoveVertices()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t n = 2 * N;
			if (k < N)
				return false;
			
			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(k);
			vertexHandler.PermutationMatrix(perm.indices());
			invG = perm.transpose() * invG * perm;
			
			matrix_t<n, n> S = invG.template bottomRightCorner<n, n>();
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), -N) * configSpace.RemovalFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor * S.determinant();
			if (acceptRatio < 0.0)
			{
				std::cout << "RemoveVertices(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "IndexBuffer:" << std::endl;
				vertexHandler.PrintIndexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				invSolver.compute(S);
				matrix_t<Eigen::Dynamic, Eigen::Dynamic> invS = invSolver.inverse();
				invG.topLeftCorner(k - n, k - n).noalias() -= invG.topRightCorner(k - n, n) * invS * invG.bottomLeftCorner(n, k - n);
				invG.conservativeResize(k - n, k - n);
				
				vertexHandler.RemoveBufferedVertices();
				return true;
			}
			else
			{
				invG = perm * invG * perm.transpose();
				return false;
			}
		}
		
		template<int_t N>
		bool RemoveVerticesWithWorms()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			uint_t l = 2 * vertexHandler.Worms();
			const uint_t n = 2 * N;
			if (k < n)
				return false;

			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(k);
			vertexHandler.PermutationMatrix(perm.indices());
			invG = perm.transpose() * invG * perm;
			wormU = perm.transpose() * wormU;
			wormV = wormV * perm;
						
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> u(k - n, n + l), v(n + l, k - n), a(n + l, n + l);
			vertexHandler.WoodburyRemoveVertices(u, v, a, perm.indices());
			u.rightCols(l) = wormU.topRows(k - n);
			v.bottomRows(l) = wormV.leftCols(k - n);
			a.bottomRightCorner(l, l) = wormA;
			a.topRightCorner(n, l) = wormU.template bottomRows<n>();
			a.bottomLeftCorner(l, n) = wormV.template rightCols<n>();
			
			matrix_t<n, n> S = invG.template bottomRightCorner<n, n>();
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> newInvG = invG.topLeftCorner(k - n, k - n);
			
			invSolver.compute(S);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> invS = invSolver.inverse();
			newInvG.noalias() -= invG.topRightCorner(k - n, n) * invS * invG.bottomLeftCorner(n, k - n);
			value_t newDetWormS = 1.0 / (wormA - wormV.leftCols(k - n) * newInvG * wormU.topRows(k - n)).determinant();
			
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), -N) * configSpace.RemovalFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor / newDetWormS / (a - v * newInvG * u).determinant();
			if (acceptRatio < 0.0)
			{
				std::cout << "RemoveVerticesWithWorm(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "IndexBuffer:" << std::endl;
				vertexHandler.PrintIndexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				invG.resize(k - n, k - n);
				invG = newInvG;
				
				wormU.conservativeResize(k - n, l);
				wormV.conservativeResize(l, k - n);
				detWormS = newDetWormS;
				
				vertexHandler.RemoveBufferedVertices();
				return true;
			}
			else
			{
				invG = perm * invG * perm.transpose();
				wormU = perm * wormU;
				wormV = wormV * perm.transpose();
				return false;
			}
		}
		
		template<int_t N>
		bool AddWorms(double preFactor)
		{
			uint_t k = 2 * vertexHandler.Vertices();
			uint_t l = 2 * vertexHandler.Worms();
			const uint_t n = 2 * N;
			if (l + n > 2 * maxWorms)
				return false;
			
			wormU.conservativeResize(k, l + n);
			wormV.conservativeResize(l + n, k);
			wormA.conservativeResize(l + n, l + n);
			vertexHandler.WoodburyAddWorm(wormU, wormV, wormA);

			matrix_t<Eigen::Dynamic, Eigen::Dynamic> invS = wormA;
			invS.noalias() -= wormV * invG * wormU;
			value_t detInvS = invS.determinant();
			value_t detRatio = detInvS * (l > 0 ? detWormS : 1.0);
			value_t acceptRatio = preFactor * detRatio * vertexHandler.VertexBufferParity();
			if (acceptRatio < 0.0)
			{
				std::cout << "AddWorm(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "VertexBuffer:" << std::endl;
				vertexHandler.PrintVertexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				detWormS = 1.0 / detInvS;
				vertexHandler.AddBufferedWorms();
				return true;
			}
			else
			{
				wormU.conservativeResize(k, l);
				wormV.conservativeResize(l, k);
				wormA.conservativeResize(l, l);
				return false;
			}
		}
		
		template<int_t N>
		bool RemoveWorms(double preFactor)
		{
			uint_t k = 2 * vertexHandler.Vertices();
			uint_t l = 2 * vertexHandler.Worms();
			const uint_t n = 2 * N;
			if (l < n)
				return false;

			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(l);
			vertexHandler.PermutationMatrix(perm.indices());

			wormU = wormU * perm;
			wormV = perm.transpose() * wormV;
			wormA = perm.transpose() * wormA * perm;

			matrix_t<Eigen::Dynamic, Eigen::Dynamic> invS(l - n, l - n);
			value_t detInvS;
			value_t detRatio;
			if (l - n > 0)
			{
				invS = wormA.topLeftCorner(l - n, l - n);
				invS.noalias() -= wormV.topRows(l - n) * invG * wormU.leftCols(l - n);
				detInvS = invS.determinant();
				detRatio = detInvS * detWormS;
			}
			else
			{
				detRatio = detWormS;
			}
			value_t acceptRatio = preFactor * detRatio * vertexHandler.WormIndexBufferParity();
			if (acceptRatio < 0.0)
			{
				std::cout << "RemoveWorm(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "IndexBuffer:" << std::endl;
				vertexHandler.PrintIndexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				if (l - n > 0)
				{
					wormU.conservativeResize(k, l - n);
					wormV.conservativeResize(l - n, k);
					wormA.conservativeResize(l - n, l - n);
					detWormS = 1.0 / detInvS;
				}

				vertexHandler.RemoveBufferedWorms();
				return true;
			}
			else
			{
				wormU = wormU * perm.transpose();
				wormV = perm * wormV;
				wormA = perm * wormA * perm.transpose();
				return false;
			}
		}
		
		bool ShiftWorm()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			uint_t l = 2 * vertexHandler.Worms();
			
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> shiftedWormU(wormU.rows(), wormU.cols());
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> shiftedWormV(wormV.rows(), wormV.cols());
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> shiftedWormA(wormA.rows(), wormA.cols());
			vertexHandler.ShiftWorm();
			vertexHandler.WoodburyWorm(shiftedWormU, shiftedWormV, shiftedWormA);
				
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> shiftedInvS = shiftedWormA - shiftedWormV * invG * shiftedWormU;
			value_t detShiftedInvS = shiftedInvS.determinant();
			value_t acceptRatio = detShiftedInvS * detWormS * vertexHandler.WormShiftParity();
			if (acceptRatio < 0.0)
			{
				std::cout << "WormShift: AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
			}
			if (configSpace.rng() < acceptRatio)
			{
				detWormS = 1.0 / detShiftedInvS;
				wormU = shiftedWormU;
				wormV = shiftedWormV;
				wormA = shiftedWormA;
				return true;
			}
			else
			{
				vertexHandler.UndoWormShift();
				return false;
			}
		}
		
		void SymmetrizeInvG()
		{
			SymmetrizeMatrix(invG);
		}

		void StabalizeInvG()
		{
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> G(invG.rows(), invG.cols());
			vertexHandler.PropagatorMatrix(G);
			invSolver.compute(G);
			invG = invSolver.inverse();
		}
		
		void StabilizeInvG(value_t& avgError, value_t& maxError)
		{
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> G(invG.rows(), invG.cols());
			vertexHandler.PropagatorMatrix(G);
			invSolver.compute(G);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> stabInvG = invSolver.inverse();
			avgError = 0.0;
			maxError = 0.0;
			value_t N = stabInvG.rows() * stabInvG.rows();
			for (uint_t i = 0; i < stabInvG.rows(); ++i)
			{
				for (uint_t j = 0; j < stabInvG.cols(); ++j)
				{
					value_t err = std::abs(invG(i, j) - stabInvG(i, j));
					if (err > maxError)
						maxError = err;
					avgError += err / N;
				}
			}
			/*
			if (avgError > std::pow(10.0, -3.0))
			{
				std::cout << "Error: " << avgError << std::endl;
				std::cout << "invG:" << std::endl;
				PrintMatrix(invG);
				std::cout << "stabInvG:" << std::endl;
				PrintMatrix(stabInvG);
				std::cin.get();
			}
			*/
			invG = stabInvG;
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
			d.write(detWormS);
			vertexHandler.Serialize(d);
		}
		
		void Serialize(idump& d)
		{
			d.read(detWormS);
			vertexHandler.Serialize(d);
			invG.resize(2 * vertexHandler.Vertices(), 2 * vertexHandler.Vertices());
			StabalizeInvG();
			if (vertexHandler.Worms() > 0)
			{
				wormU.resize(2 * vertexHandler.Vertices(), 2 * vertexHandler.Worms());
				wormV.resize(2 * vertexHandler.Worms(), 2 * vertexHandler.Vertices());
				wormA.resize(2 * vertexHandler.Worms(), 2 * vertexHandler.Worms());
				vertexHandler.WoodburyWorm(wormU, wormV, wormA);
			}
		}
		
	private:
		ConfigSpace_t& configSpace;
		VertexHandler_t vertexHandler;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> invG;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormU;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormV;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormA;
		value_t detWormS;
		Eigen::FullPivLU< matrix_t<Eigen::Dynamic, Eigen::Dynamic> > invSolver;
		uint_t maxWorms = 2;
};