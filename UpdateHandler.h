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
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice->Bonds(), N) * configSpace.AdditionFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor * invS.determinant();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "AddVertices(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "VertexBuffer:" << std::endl;
				vertexHandler.PrintVertexBuffer();
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
				
				vertexHandler.AddBufferedVertices();
				return true;
			}
			return false;
		}
		
		template<int_t N, int_t W>
		bool AddVerticesWithWorms()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t n = 2 * N;
			const uint_t l = 2 * W;
			
			matrix_t<Eigen::Dynamic, l + n> u(k, l + n);
			matrix_t<l + n, Eigen::Dynamic> v(l + n, k);
			matrix_t<l + n, l + n> a(l + n, l + n);
			vertexHandler.WoodburyAddVertices(u, v, a);
			u.rightCols(l) = wormU;
			v.bottomRows(l) = wormV;
			a.template bottomRightCorner<l, l>() = wormA;

			matrix_t<l + n, l + n> invS = a;
			invS.noalias() -= v * invG * u;
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice->Bonds(), N) * configSpace.AdditionFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor * invS.determinant() * detWormS;
			if (print && acceptRatio < 0.0)
			{
				std::cout << "AddVerticesWithWorm(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "Worms:" << std::endl;
				vertexHandler.PrintWormVertices();
				std::cout << "VertexBuffer:" << std::endl;
				vertexHandler.PrintVertexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				wormU.conservativeResize(k + n, l);
				wormU.bottomRows(n) = a.template topRightCorner<n, l>();
				wormV.conservativeResize(l, k + n);
				wormV.rightCols(n) = a.template bottomLeftCorner<l, n>();
				
				matrix_t<n, n> invS = a.template topLeftCorner<n, n>();
				matrix_t<Eigen::Dynamic, n> invGu = invG * u.topLeftCorner(k, n);
				invS.noalias() -= v.topLeftCorner(n, k) * invGu;
				
				matrix_t<n, n> S = invS.inverse();
				matrix_t<n, Eigen::Dynamic> R = -S * v.topLeftCorner(n, k) * invG;
				
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
			if (k < n)
				return false;
			
			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(k);
			vertexHandler.PermutationMatrix(perm.indices());
			invG = perm.transpose() * invG * perm;
			
			matrix_t<n, n> S = invG.template bottomRightCorner<n, n>();
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice->Bonds(), -N) * configSpace.RemovalFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor * S.determinant();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "RemoveVertices(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "IndexBuffer:" << std::endl;
				vertexHandler.PrintIndexBuffer();
			}
			if (configSpace.rng() < acceptRatio)
			{
				matrix_t<n, n> invS = S.inverse();
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
		
		template<int_t N, int_t W>
		bool RemoveVerticesWithWorms()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t n = 2 * N;
			const uint_t l = 2 * W;
			if (k < n)
				return false;

			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(k);
			vertexHandler.PermutationMatrix(perm.indices());
			invG = perm.transpose() * invG * perm;
			wormU = perm.transpose() * wormU;
			wormV = wormV * perm;
						
			matrix_t<Eigen::Dynamic, n + l> u(k - n, n + l);
			matrix_t<n + l, Eigen::Dynamic> v(n + l, k - n);
			matrix_t<n + l, n + l> a(n + l, n + l);
			vertexHandler.WoodburyRemoveVertices(u, v, a, perm.indices());
			u.rightCols(l) = wormU.topRows(k - n);
			v.bottomRows(l) = wormV.leftCols(k - n);
			a.template bottomRightCorner<l, l>() = wormA;
			a.template topRightCorner<n, l>() = wormU.template bottomRows<n>();
			a.template bottomLeftCorner<l, n>() = wormV.template rightCols<n>();
			
			matrix_t<n, n> S = invG.template bottomRightCorner<n, n>();
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> newInvG = invG.topLeftCorner(k - n, k - n);
			
			newInvG.noalias() -= invG.topRightCorner(k - n, n) * S.inverse() * invG.bottomLeftCorner(n, k - n);
			value_t newDetWormS = 1.0 / (wormA - wormV.leftCols(k - n) * newInvG * wormU.topRows(k - n)).determinant();
			
			value_t preFactor = std::pow(-configSpace.beta * configSpace.V * configSpace.lattice->Bonds(), -N) * configSpace.RemovalFactorialRatio(k / 2, N);
			value_t acceptRatio = preFactor / newDetWormS / (a - v * newInvG * u).determinant();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "RemoveVerticesWithWorm(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "Worms:" << std::endl;
				vertexHandler.PrintWormVertices();
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
		
		template<int_t N, int_t W>
		bool AddWorms(double preFactor)
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t n = 2 * N;
			const uint_t l = 2 * W;
			if (l + n > 2 * maxWorms)
				return false;
			
			wormU.conservativeResize(k, l + n);
			wormV.conservativeResize(l + n, k);
			wormA.conservativeResize(l + n, l + n);
			vertexHandler.WoodburyAddWorm(wormU, wormV, wormA);

			matrix_t<l + n, l + n> invS = wormA;
			invS.noalias() -= wormV * invG * wormU;
			value_t detInvS = invS.determinant();
			value_t detRatio = detInvS * (l > 0 ? detWormS : 1.0);
			value_t acceptRatio = preFactor * detRatio * vertexHandler.VertexBufferParity();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "AddWorm(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				matrix_t<Eigen::Dynamic, Eigen::Dynamic> G(invG.rows(), invG.cols());
				vertexHandler.PropagatorMatrix(G);
				std::cout << MatrixCondition(G) << std::endl;
				/*
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "Worms:" << std::endl;
				vertexHandler.PrintWormVertices();
				std::cout << "VertexBuffer:" << std::endl;
				vertexHandler.PrintVertexBuffer();
				*/
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
		
		template<int_t N, int_t W>
		bool RemoveWorms(double preFactor)
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t n = 2 * N;
			const uint_t l = 2 * W;
			if (l < n)
				return false;

			Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(l);
			vertexHandler.PermutationMatrix(perm.indices());

			wormU = wormU * perm;
			wormV = perm.transpose() * wormV;
			wormA = perm.transpose() * wormA * perm;

			matrix_t<l - n, l - n> invS(l - n, l - n);
			value_t detInvS;
			value_t detRatio;
			if (l - n > 0)
			{
				invS = wormA.template topLeftCorner<l - n, l - n>();
				invS.noalias() -= wormV.topRows(l - n) * invG * wormU.leftCols(l - n);
				detInvS = invS.determinant();
				detRatio = detInvS * detWormS;
			}
			else
			{
				detRatio = detWormS;
			}
			value_t acceptRatio = preFactor * detRatio * vertexHandler.WormIndexBufferParity();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "RemoveWorm(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "Worms:" << std::endl;
				vertexHandler.PrintWormVertices();
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
		
		template<int_t W>
		bool ShiftWorm()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t l = 2 * W;
			
			matrix_t<Eigen::Dynamic, l> shiftedWormU(wormU.rows(), wormU.cols());
			matrix_t<l, Eigen::Dynamic> shiftedWormV(wormV.rows(), wormV.cols());
			matrix_t<l, l> shiftedWormA(wormA.rows(), wormA.cols());
			vertexHandler.ShiftWorm();
			vertexHandler.WoodburyWorm(shiftedWormU, shiftedWormV, shiftedWormA);
				
			matrix_t<l, l> shiftedInvS = shiftedWormA;
			shiftedInvS.noalias() -= shiftedWormV * invG * shiftedWormU;
			value_t detShiftedInvS = shiftedInvS.determinant();
			value_t acceptRatio = detShiftedInvS * detWormS * vertexHandler.WormShiftParity();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "WormShift: AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "Vertices:" << std::endl;
				vertexHandler.PrintVertices();
				std::cout << "Worms:" << std::endl;
				vertexHandler.PrintWormVertices();
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
		
		template<typename Matrix>
		value_t MatrixCondition(Matrix& M)
		{
			if (M.rows() == 0)
				return 0.0;
			Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
			//std::cout << svd.singularValues()(M.rows()-1) << std::endl;
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

		double PropagatorMatrixCondition()
		{
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> G(2 * vertexHandler.Vertices() + 2 * vertexHandler.Worms(), 2 * vertexHandler.Vertices() + 2 * vertexHandler.Worms());
			vertexHandler.FullPropagatorMatrix(G);
			return MatrixCondition(G);
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
			StabilizeInvG();
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
		Matrices lastCheckpoint;
		VertexHandler_t vertexHandler;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> invG;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormU;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormV;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> wormA;
		value_t detWormS;
		uint_t maxWorms = 2;
		bool print = false;
};