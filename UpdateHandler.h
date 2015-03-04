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

enum UpdateFlag {NormalUpdate, NoUpdate, ForceUpdate};

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
			value_t det;
			return AddVertices<N>(preFactor, isWorm, det, UpdateFlag::NormalUpdate);
		}

		template<int_t N>
		bool AddVertices(value_t preFactor, bool isWorm, value_t& det, UpdateFlag flag)
		{
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;
			GeMatrix u(k, n);
			GeMatrix v(n, k);
			GeMatrix a(n, n);
			vertexHandler.WoodburyAddVertices(u, v, a);

			GeMatrix invGu = invG * u;
			GeMatrix S = a - v * invGu;
			IndexVector pivS(n);

			MatrixOperation<value_t, 2*N> matop;
			value_t acceptRatio;
			if (flag == UpdateFlag::NormalUpdate)
			{
				det = matop.Determinant(S, pivS);
				acceptRatio = preFactor * det;
			}
			else if (flag == UpdateFlag::NoUpdate)
			{
				det = matop.Determinant(S, pivS);
				acceptRatio = 0.0;
			}
			else if (flag == UpdateFlag::ForceUpdate)
			{
				acceptRatio = 1.0;
			}
			if (print && acceptRatio < 0.0)
			{
				std::cout << "AddVertices(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
			}
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
				vertexHandler.AddBufferedVertices(isWorm);
			}
			return false;
		}

		template<int_t N>
		bool RemoveVertices(value_t preFactor, bool isWorm)
		{
			value_t det;
			return RemoveVertices<N>(preFactor, isWorm, det, UpdateFlag::NormalUpdate);
		}

		template<int_t N>
		bool RemoveVertices(value_t preFactor, bool isWorm, value_t& det, UpdateFlag flag)
		{
			if (isWorm && vertexHandler.Worms() < N)
				return false;
			if ((!isWorm) && vertexHandler.Vertices() < N)
				return false;
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;

			GeMatrix perm(k, k);
			vertexHandler.PermutationMatrix(perm, isWorm);
			GeMatrix invGp = invG * perm;
			GeMatrix pInvGp = flens::transpose(perm) * invGp;
			
			const flens::Underscore<IndexType> _;
			typename GeMatrix::View S = pInvGp(_(k - n + 1, k), _(k - n + 1, k));
			IndexVector pivS(n);
			MatrixOperation<value_t, 2*N> matop;
			value_t acceptRatio;
			if (flag == UpdateFlag::NormalUpdate)
			{
				det = matop.Determinant(S, pivS);
				acceptRatio = preFactor * det;
			}
			else if (flag == UpdateFlag::NoUpdate)
			{
				det = matop.Determinant(S, pivS);
				acceptRatio = 0.0;
			}
			else if (flag == UpdateFlag::ForceUpdate)
			{
				acceptRatio = 1.0;
			}
			if (print && acceptRatio < 0.0)
			{
				std::cout << "RemoveVertices(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
			}
			if (configSpace.rng() < acceptRatio)
			{
				typename GeMatrix::View P = pInvGp(_(1, k - n), _(1, k - n));
				typename GeMatrix::View Q = pInvGp(_(1, k - n), _(k - n + 1, k));
				typename GeMatrix::View R = pInvGp(_(k - n + 1, k), _(1, k - n));
				matop.Inverse(S, pivS);
				invG.resize(k - n, k - n);
				GeMatrix SR = S * R;
				invG = P - Q * SR;
				
				vertexHandler.RemoveBufferedVertices(isWorm);
				return true;
			}
			else
			{
				return false;
			}
		}

		bool OpenUpdate(value_t preFactor)
		{
			value_t acceptRatio = preFactor;
			if (configSpace.rng() < acceptRatio)
			{
				if (vertexHandler.Vertices() == 0)
					return false;
				else
				{
					vertexHandler.OpenUpdate();
					return true;
				}
			}
			else
			{
				return false;
			}
		}

		bool CloseUpdate(value_t preFactor)
		{
			value_t acceptRatio = preFactor;
			if (configSpace.rng() < acceptRatio)
			{
				vertexHandler.CloseUpdate();
				return true;
			}
			else
			{
				return false;
			}
		}
		
	/*	
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

			value_t preFactorRem, preFactorAdd;
			value_t m = configSpace.lattice->Sites();
			if (W == 1)
			{
				preFactorRem = 1.0 / (configSpace.lattice->Sites() * m * configSpace.beta * configSpace.zeta2);
				preFactorAdd = 1.0 / preFactorRem;
			}
			else if (W == 2)
			{
				preFactorRem = 1.0 / (configSpace.lattice->Sites() * m * m * m * configSpace.beta * configSpace.zeta4);
				preFactorAdd = 1.0 / preFactorRem;
			}

			value_t acceptRatio = detShiftedInvS / detInvS * vertexHandler.WormShiftParity();
			//value_t acceptRatio = std::min({std::abs(preFactorRem * detShiftedInvS), 1.0}) * std::min({std::abs(preFactorAdd / detInvS), 1.0});
			if (print && acceptRatio < 0.0)
			{
				std::cout << "WormShift: AcceptRatio: " << acceptRatio << std::endl;
			}
			if (configSpace.rng() < acceptRatio)
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
		*/
		
		template<int_t W>
		bool ShiftWorm()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			const uint_t l = 2 * W;

			vertexHandler.ShiftWormToBuffer();
			value_t preFactorRem, preFactorAdd;
			value_t m = configSpace.lattice->Sites();
			if (W == 1)
			{
				preFactorRem = 1.0 / (configSpace.lattice->Sites() * m * configSpace.beta * configSpace.zeta2);
				preFactorAdd = 1.0 / preFactorRem;
			}
			else if (W == 2)
			{
				preFactorRem = 1.0 / (configSpace.lattice->Sites() * m * m * m * configSpace.beta * configSpace.zeta4);
				preFactorAdd = 1.0 / preFactorRem;
			}
			if (RemoveVertices<W>(preFactorRem * vertexHandler.WormIndexBufferParity(), true))
			{
				if (AddVertices<W>(preFactorAdd * vertexHandler.VertexBufferParity(), true))
				{
					return true;
				}
				else
				{
					vertexHandler.template RestoreAfterShift<W>();
					value_t det;
					AddVertices<W>(preFactorAdd * vertexHandler.VertexBufferParity(), true, det, UpdateFlag::ForceUpdate);
					return false;
				}
			}
			return false;
		}
		

		void Clear()
		{
			invG.resize(0, 0);
			vertexHandler.Clear();
		}
		
		template<typename Matrix>
		value_t MatrixCondition(Matrix& M)
		{
			if (M.rows() == 0)
				return 0.0;
			Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
			return svd.singularValues()(0) / svd.singularValues()(M.rows()-1);
		}

		void SymmetrizeInvG()
		{
			SymmetrizeMatrix(invG);
		}

		value_t StabilizeInvG()
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
		
		value_t StabilizeInvG(value_t& avgError)
		{
			if (vertexHandler.Vertices() == 0)
				return 0.0;
			GeMatrix stabInvG(invG.numRows(), invG.numCols());
			vertexHandler.PropagatorMatrix(stabInvG);
			MatrixOperation<value_t, 0> matop;
			matop.Inverse(stabInvG);
			avgError = 0.0;
			value_t N = stabInvG.numRows() * stabInvG.numRows();
			for (uint_t i = 1; i <= stabInvG.numRows(); ++i)
			{
				for (uint_t j = 1; j <= stabInvG.numCols(); ++j)
				{
					value_t err = std::abs(invG(i, j) - stabInvG(i, j));
					avgError += err / N;
				}
			}
			invG = stabInvG;
			return 0.0;
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
			invG.resize(2 * (vertexHandler.Vertices() + vertexHandler.Worms()), 2 * (vertexHandler.Vertices() + vertexHandler.Worms()));
			StabilizeInvG();
		}
		
	private:
		ConfigSpace_t& configSpace;
		VertexHandler_t vertexHandler;
		GeMatrix invG;
		uint_t maxWorms = 2;
		bool print = true;
};
