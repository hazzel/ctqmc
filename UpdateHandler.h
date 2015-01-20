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
#include "MatrixOperation.h"

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
		typedef flens::DenseVector< flens::Array<value_t> > GeVector;
		
		struct Matrices
		{
			GeMatrix invG;
			GeMatrix wormU;
			GeMatrix wormV;
			GeMatrix wormA;
			value_t detWormS;
		};
		
		UpdateHandler(ConfigSpace_t& configSpace)
			: configSpace(configSpace), vertexHandler(VertexHandler_t(configSpace)), invG(GeMatrix(0, 0))
		{}
		
		void Init()
		{
			for (int_t n = 1; n <= 2; ++n)
				prefactor.push_back(std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), n));
		}
		
		//Computes SVD
		template<typename Matrix, typename Vector>
		void SVD(Matrix& M, Vector& S)
		{
			GeMatrix U(M.numRows(), M.numCols());
			GeMatrix VT(M.numRows(), M.numCols());
			flens::lapack::svd(flens::lapack::SVD::Job::All, flens::lapack::SVD::Job::All, M, S, U, VT);
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
			
			MatrixOperation<value_t, 2*N> matop;
			value_t detS = matop.Determinant(S, pivS);
			value_t preFactor = configSpace.AdditionFactorialRatio(k / 2, N) * prefactor[N-1];
			value_t acceptRatio = preFactor * detS;
			if (acceptRatio < 0.0)
				std::cout << "AddVertices: AcceptRatio: " << acceptRatio << std::endl;
			if (configSpace.rng() < acceptRatio)
			{
				matop.Inverse(S, pivS);
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
		bool AddVerticesWithWorms()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			uint_t l = 2 * vertexHandler.Worms();
			const uint_t n = 2 * N;
		
			const flens::Underscore<IndexType> _;
			GeMatrix u(k, n + l), v(n + l, k), a(n + l, n + l);
			vertexHandler.WoodburyAddVertices(u, v, a);
			u(_, _(n + 1, n + l)) = wormU;
			v(_(n + 1, n + l), _) = wormV;
			a(_(n + 1, n + l), _(n + 1, n + l)) = wormA;

			GeMatrix invGu = invG * u;
			GeMatrix invS = a - v * invGu;
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), N) * configSpace.AdditionFactorialRatio(k / 2, N);
			MatrixOperation<value_t, 0> matop;
			value_t acceptRatio = preFactor * matop.Determinant(invS) * detWormS;
			if (acceptRatio < 0.0)
			{
				std::cout << "AddVerticesWithWorm(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				vertexHandler.PrintVertices();
				vertexHandler.PrintVertexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				wormU.resize(k + n, l);
				wormU(_(1, k), _) = u(_, _(n + 1, n + l));
				wormU(_(k + 1, k + n), _) = a(_(1, n), _(n + 1, n + l));
				wormV.resize(l, k + n);
				wormV(_, _(1, k)) = v(_(n + 1, n + l), _);
				wormV(_, _(k + 1, k + n)) = a(_(n + 1, n + l), _(1, n));

				typename GeMatrix::View u_view = u(_(1, k), _(1, n));
				typename GeMatrix::View v_view = v(_(1, n), _(1, k));
				typename GeMatrix::View a_view = a(_(1, n), _(1, n));

				GeMatrix invGu = invG * u_view;
				GeMatrix vinvG = v_view * invG;
				GeMatrix S = a_view - v_view * invGu;
				
				matop.Inverse(S);
				GeMatrix R = -S * vinvG;
				GeMatrix P = invG - invGu * R;
				
				invG.resize(k + n, k + n);
				invG(_(1, k), _(1, k)) = P;
				invG(_(1, k), _(k + 1, k + n)) = -invGu * S;
				invG(_(k + 1, k + n), _(1, k)) = R;
				invG(_(k + 1, k + n), _(k + 1, k + n)) = S;
				
				GeMatrix invGwU = invG * wormU;
				GeMatrix wormS = wormA - wormV * invGwU;

				detWormS = 1.0 / matop.Determinant(wormS);
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
			if (k < n)
				return false;
		
			GeMatrix perm(k, k);
			vertexHandler.PermutationMatrix(perm);
			GeMatrix invGp = invG * perm;
			GeMatrix pInvGp = flens::transpose(perm) * invGp;
			
			const flens::Underscore<IndexType> _;
			typename GeMatrix::View S = pInvGp(_(k - n + 1, k), _(k - n + 1, k));
			IndexVector pivS(n);
			MatrixOperation<value_t, 2*N> matop;
			value_t detS = matop.Determinant(S, pivS);
			value_t preFactor = configSpace.RemovalFactorialRatio(k / 2, N) / prefactor[N-1];
			value_t acceptRatio = preFactor * detS;
			if (acceptRatio < 0.0)
				std::cout << "RemoveVertex: AcceptRatio" << acceptRatio << std::endl;
			if (configSpace.rng() < acceptRatio)
			{
				typename GeMatrix::View P = pInvGp(_(1, k - n), _(1, k - n));
				typename GeMatrix::View Q = pInvGp(_(1, k - n), _(k - n + 1, k));
				typename GeMatrix::View R = pInvGp(_(k - n + 1, k), _(1, k - n));
				matop.Inverse(S, pivS);
				invG.resize(k - n, k - n);
				GeMatrix SR = S * R;
				invG = P - Q * SR;
				
				vertexHandler.RemoveBufferedVertices();
				return true;
			}
			else
			{
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
			
			const flens::Underscore<IndexType> _;
			GeMatrix perm(k, k);
			std::vector<uint_t> perm_indices(k);
			vertexHandler.PermutationMatrix(perm, perm_indices);
			GeMatrix invGp = invG * perm;
			GeMatrix pInvGp = flens::transpose(perm) * invGp;
			GeMatrix wormUp = flens::transpose(perm) * wormU;
			typename GeMatrix::View wU = wormUp(_(1, k - n), _);
			GeMatrix wormVp = wormV * perm;
			typename GeMatrix::View wV = wormVp(_, _(1, k - n));
			
			GeMatrix u(k - n, n + l), v(n + l, k - n), a(n + l, n + l);
			vertexHandler.WoodburyRemoveVertices(u, v, a, perm_indices);
			u(_, _(n + 1, n + l)) = wU;
			v(_(n + 1, n + l), _) = wV;
			a(_(n + 1, n + l), _(n + 1, n + l)) = wormA;
			a(_(1, n), _(n + 1, n + l)) = wormUp(_(k - n + 1, k), _);
			a(_(n + 1, n + l), _(1, n)) = wormVp(_, _(k - n + 1, k));
			
			typename GeMatrix::View P = pInvGp(_(1, k - n), _(1, k - n));
			typename GeMatrix::View Q = pInvGp(_(1, k - n), _(k - n + 1, k));
			typename GeMatrix::View R = pInvGp(_(k - n + 1, k), _(1, k - n));
			typename GeMatrix::View S = pInvGp(_(k - n + 1, k), _(k - n + 1, k));
			MatrixOperation<value_t, 0> matop;
				
			matop.Inverse(S);
			GeMatrix SR = S * R;
			GeMatrix newInvG = P - Q * SR;
			
			GeMatrix newInvGwU = newInvG * wU;
			GeMatrix invWormS = wormA - wV * newInvGwU;
			value_t newDetWormS = 1.0 / matop.Determinant(invWormS);
			
			GeMatrix newInvGu = newInvG * u;
			GeMatrix invS = a - v * newInvGu;
			
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice.Bonds(), -N) * configSpace.RemovalFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor / newDetWormS / matop.Determinant(invS);;
			if (acceptRatio < 0.0)
			{
				std::cout << "RemoveVerticesWithWorm(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "IndexBuffer:" << std::endl;
				vertexHandler.PrintIndexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				invG.resize(k - n, k - n);
				invG = newInvG;
				
				wormU.resize(k - n, l);
				wormU = wU;
				wormV.resize(l, k - n);
				wormV = wV;
				detWormS = newDetWormS;
				
				vertexHandler.RemoveBufferedVertices();
				return true;
			}
			else
			{
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
			
			const flens::Underscore<IndexType> _;
			GeMatrix newWormU(k, l + n);
			GeMatrix newWormV(l + n, k);
			GeMatrix newWormA(l + n, l + n);
			if (l > 0)
			{
				newWormU(_, _(1, l)) = wormU;
				newWormV(_(1, l), _) = wormV;
				newWormA(_(1, l), _(1, l)) = wormA;
			}
			vertexHandler.WoodburyAddWorm(newWormU, newWormV, newWormA);

			MatrixOperation<value_t, 0> matop;
			GeMatrix invGwU = invG * newWormU;
			GeMatrix invS = newWormA - newWormV * invGwU;
			value_t detInvS = matop.Determinant(invS);
			value_t detRatio = detInvS * (l > 0 ? detWormS : 1.0);
			value_t acceptRatio = preFactor * detRatio * vertexHandler.VertexBufferParity();
			if (acceptRatio < 0.0)
			{
				std::cout << "AddWorm(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				vertexHandler.PrintVertices();
				vertexHandler.PrintVertexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				detWormS = 1.0 / detInvS;
				
				wormU.resize(k, l + n);
				wormU = newWormU;
				wormV.resize(l + n, k);
				wormV = newWormV;
				wormA.resize(l + n, l + n);
				wormA = newWormA;
				
				vertexHandler.AddBufferedWorms();
				return true;
			}
			else
			{
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

			const flens::Underscore<IndexType> _;
			GeMatrix perm(l, l);
			vertexHandler.PermutationMatrix(perm);

			GeMatrix wormUp = wormU * perm;
			GeMatrix wormVp = flens::transpose(perm) * wormV;
			GeMatrix O = wormA * perm;
			GeMatrix wormAp = flens::transpose(perm) * O;

			value_t detInvS;
			value_t detRatio;
			typename GeMatrix::View newWormU = wormUp(_, _(1, l - n));
			typename GeMatrix::View newWormV = wormVp(_(1, l - n), _);
			typename GeMatrix::View newWormA = wormAp(_(1, l - n), _(1, l - n));
			MatrixOperation<value_t, 0> matop;
			if (l - n > 0)
			{
				GeMatrix invGwU = invG * newWormU;
				GeMatrix invS = newWormA - newWormV * invGwU;
				detInvS = matop.Determinant(invS);
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
				vertexHandler.PrintVertices();
				std::cout << "IndexBuffer:" << std::endl;
				vertexHandler.PrintIndexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				if (l - n > 0)
				{
					wormU.resize(k, l - n);
					wormU = newWormU;
					wormV.resize(l - n, k);
					wormV = newWormV;
					wormA.resize(l - n, l - n);
					wormA = newWormA;
					detWormS = 1.0 / detInvS;
				}

				vertexHandler.RemoveBufferedWorms();
				return true;
			}
			else
			{
				return false;
			}
		}
		
		bool ShiftWorm()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			uint_t l = 2 * vertexHandler.Worms();
			
			GeMatrix shiftedWormU(wormU.numRows(), wormU.numCols());
			GeMatrix shiftedWormV(wormV.numRows(), wormV.numCols());
			GeMatrix shiftedWormA(wormA.numRows(), wormA.numCols());
			vertexHandler.ShiftWorm();
			vertexHandler.WoodburyWorm(shiftedWormU, shiftedWormV, shiftedWormA);
			
			MatrixOperation<value_t, 0> matop;
			GeMatrix invGwU = invG * shiftedWormU;
			GeMatrix shiftedInvS = shiftedWormA - shiftedWormV * invGwU;
			value_t detShiftedInvS = matop.Determinant(shiftedInvS);
			value_t acceptRatio = detShiftedInvS * detWormS * vertexHandler.WormShiftParity();
			if (acceptRatio < 0.0)
			{
				std::cout << "WormShift: AcceptRatio: " << acceptRatio << std::endl;
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

		value_t StabalizeInvG()
		{
			if (vertexHandler.Vertices() == 0)
				return 0.0;
			GeMatrix stabInvG(invG.numRows(), invG.numCols());
			vertexHandler.PropagatorMatrix(stabInvG);
			MatrixOperation<value_t, 0> matop;
			matop.Inverse(stabInvG);
			invG = stabInvG;
			return 0.0;
		}
		
		value_t StabilizeInvG(value_t& avgError, value_t& relError)
		{
			if (vertexHandler.Vertices() == 0)
				return 0.0;
			GeMatrix stabInvG(invG.numRows(), invG.numCols());
			vertexHandler.PropagatorMatrix(stabInvG);
			MatrixOperation<value_t, 0> matop;
			matop.Inverse(stabInvG);
			avgError = 0.0;
			relError = 0.0;
			value_t N = stabInvG.numRows() * stabInvG.numRows();
			for (uint_t i = 1; i <= stabInvG.numRows(); ++i)
			{
				for (uint_t j = 1; j <= stabInvG.numCols(); ++j)
				{
					value_t err = std::abs(invG(i, j) - stabInvG(i, j));
					avgError += err / N;
					relError += std::abs(stabInvG(i, j));
				}
			}
			relError = avgError * N / relError;
			invG = stabInvG;
			return 0.0;
		}
		
		template<typename Matrix>
		void SymmetrizeMatrix(Matrix& M)
		{
			for (uint_t i = 1; i < M.numRows(); ++i)
			{
				for (uint_t j = 1; j < i; ++j)
				{
					value_t mean = (std::abs(M(i, j)) + std::abs(M(j, i))) / 2.0;
					M(i, j) = sgn(M(i, j)) * mean;
					M(j, i) = sgn(M(j, i)) * mean;
				}
				M(i, i) = 0.;
			}
		}
		
		void PrintPropagatorMatrix()
		{
			GeMatrix G(invG.numRows(), invG.numCols());
			vertexHandler.PropagatorMatrix(G);
			std::cout << "G:" << std::endl;
			std::cout << G << std::endl;
			std::cout << "SVD (G):" << std::endl;
			GeVector S(G.numRows());
			SVD(G, S);
			std::cout << S << std::endl;
			std::cout << "invG:" << std::endl;
			std::cout << invG << std::endl;
			std::cout << "SVD (invG):" << std::endl;
			SVD(invG, S);
			std::cout << S << std::endl;
		}
		
		void SaveCheckpoint()
		{
			lastCheckpoint.invG = invG;
			lastCheckpoint.wormU = wormU;
			lastCheckpoint.wormV = wormV;
			lastCheckpoint.wormA = wormA;
			lastCheckpoint.detWormS = detWormS;
			vertexHandler.SaveCheckpoint();
		}
		
		void RestoreCheckpoint()
		{
			invG = lastCheckpoint.invG;
			wormU = lastCheckpoint.wormU;
			wormV = lastCheckpoint.wormV;
			wormA = lastCheckpoint.wormA;
			detWormS = lastCheckpoint.detWormS;
			vertexHandler.RestoreCheckpoint();
			std::cout << std::endl;
			PrintMatrix(invG);
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
		
		GeMatrix& InvG()
		{
			return invG;
		}
		
	private:
		ConfigSpace_t& configSpace;
		Matrices lastCheckpoint;
		VertexHandler_t vertexHandler;
		GeMatrix invG;
		GeMatrix wormU;
		GeMatrix wormV;
		GeMatrix wormA;
		value_t detWormS;
		uint_t maxWorms = 2;
		std::vector<value_t> prefactor;
};