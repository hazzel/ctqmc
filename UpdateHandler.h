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
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t n = 2 * N;
			GeMatrix u(k, n);
			GeMatrix v(n, k);
			GeMatrix a(n, n);
			vertexHandler.WoodburyAddVertices(u, v, a);

			GeMatrix invGu = invG * u;
			GeMatrix S = a - v * invGu;
			IndexVector pivS(n);
			
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
				GeMatrix vinvG = v * invG;
				GeMatrix R = -S * vinvG;
				GeMatrix P = invG - invGu * R;
				
				invG.resize(k + n, k + n);
				const flens::Underscore<IndexType> _;
				invG(_(1, k), _(1, k)) = P;
				invG(_(1, k), _(k + 1, k + n)) = -invGu * S;
				invG(_(k + 1, k + n), _(1, k)) = R;
				invG(_(k + 1, k + n), _(k + 1, k + n)) = S;
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
			
			GeMatrix perm(k, k);
			vertexHandler.PermutationMatrix(perm);
			GeMatrix P = invG * perm;
			invG = flens::transpose(perm) * P;
			
			const flens::Underscore<IndexType> _;
			GeMatrix S = invG(_(k - n + 1, k), _(k - n + 1, k));
			IndexVector pivS(n);
			value_t detS = (N == 1 ? Determinant2x2(S) : Determinant(S, pivS));
			value_t preFactor = configSpace.RemovalFactorialRatio(k, N) / prefactor[N-1];
			value_t acceptRatio = preFactor * detS;
			if (acceptRatio < 0.0)
				std::cout << "RemoveVertex: AcceptRatio" << acceptRatio << std::endl;
			if (configSpace.rng() < acceptRatio)
			{
				GeMatrix P = invG(_(1, k - n), _(1, k - n));
				GeMatrix Q = invG(_(1, k - n), _(k - n + 1, k));
				GeMatrix R = invG(_(k - n + 1, k), _(1, k - n));
				if(N == 1)
					Inverse2x2(S, detS);
				else
					Inverse(S, pivS);
				
				invG.resize(k - n, k - n);
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
		
		template<int_t N>
		bool RemoveVerticesWithWorms()
		{
			return false;
			/*
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
			u.topRightCorner(k - n, l) = wormU.topRows(k - n);
			v.bottomLeftCorner(l, k - n) = wormV.leftCols(k - n);
			a.bottomRightCorner(l, l) = wormA;
			a.topRightCorner(n, l) = wormU.bottomRows(n);
			a.bottomLeftCorner(l, n) = wormV.rightCols(n);
			
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> S = invG.bottomRightCorner(n, n);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> newInvG = invG.topLeftCorner(k - n, k - n);
			newInvG.noalias() -= invG.topRightCorner(k - n, n) * S.inverse() * invG.bottomLeftCorner(n, k - n);
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
			*/
		}
		
		template<int_t N>
		bool AddWorms(double preFactor)
		{
			return false;
			/*
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
			*/
		}
		
		template<int_t N>
		bool RemoveWorms(double preFactor)
		{
			return false;
			/*
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
			*/
		}
		
		bool ShiftWorm()
		{
			return false;
			/*
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
			*/
		}
		
		void SymmetrizeInvG()
		{
			SymmetrizeMatrix(invG);
		}

		void StabalizeInvG()
		{
			/*
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> G(invG.rows(), invG.cols());
			vertexHandler.PropagatorMatrix(G);
			invSolver.compute(G);
			invG = invSolver.inverse();
			*/
		}
		
		void StabilizeInvG(value_t& avgError, value_t& maxError)
		{
			/*
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
			invG = stabInvG;
			*/
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
		GeMatrix invG;
		std::vector<value_t> prefactor;
};