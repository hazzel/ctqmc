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
#include <armadillo>
#include "measurements.h"
#include "random.h"
#include "parser.h"
#include "types.h"
#include "ConfigSpace.h"
#include "VertexHandler.h"

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
		using matrix_t = arma::mat;
		
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
			matrix_t u(k, n), v(n, k), a(n, n);
			vertexHandler.WoodburyAddVertices(u, v, a);

			matrix_t invGu = invG * u;
			matrix_t invS = a - v * invGu;

			value_t acceptRatio;
			if (flag == UpdateFlag::NormalUpdate)
			{
				det = arma::det(invS);
				acceptRatio = preFactor * det;
			}
			else if (flag == UpdateFlag::NoUpdate)
			{
				det = arma::det(invS);
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
				matrix_t S = arma::inv(invS);
				matrix_t vinvG = v * invG;
				
				matrix_t newInvG(k + n, k + n);
				matrix_t s = -S * vinvG;
				if (k != 0)
				{
					newInvG.submat(k, 0, n + k - 1, k - 1) = -S * vinvG;
					newInvG.submat(0, 0, k - 1, k - 1) = invG - invGu * newInvG.submat(k, 0, k + n - 1, k - 1);
					newInvG.submat(0, k, k - 1, k + n - 1) = -invGu * S;
				}
				newInvG.submat(k, k, k + n - 1, k + n - 1) = S;
				invG = newInvG;
				
				vertexHandler.AddBufferedVertices(isWorm);
				return true;
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

			arma::vec perm(k);
			vertexHandler.PermutationMatrix(perm, isWorm);
			matrix_t invGp(k, k);
			for (uint_t i = 0; i < k; ++i)
				for (uint_t j = 0; j < k; ++j)
					invGp(i, j) = invG(perm(i), perm(j));
			
			matrix_t S = invGp.submat(k - n, k - n, k - 1, k - 1);
			value_t acceptRatio;
			if (flag == UpdateFlag::NormalUpdate)
			{
				det = arma::det(S);
				acceptRatio = preFactor * det;
			}
			else if (flag == UpdateFlag::NoUpdate)
			{
				det = arma::det(S);
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
				matrix_t t = arma::inv(S) * invGp.submat(k - n, 0, k - 1, k - n - 1);
				invG = invGp.submat(0, 0, k - n - 1, k - n - 1) - invGp.submat(0, k - n, k - n - 1, k - 1) * t;

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
			return 0.0;
		}

		value_t StabilizeInvG()
		{
			if (invG.n_rows == 0)
				return 0.0;
			matrix_t G(invG.n_rows, invG.n_cols);
			vertexHandler.PropagatorMatrix(G);
			invG = arma::inv(G);
			return 0.0;
		}
		
		value_t StabilizeInvG(value_t& avgError)
		{
			if (invG.n_rows == 0)
				return 0.0;
			matrix_t G(invG.n_rows, invG.n_cols);
			vertexHandler.PropagatorMatrix(G);
			matrix_t stabInvG = arma::inv(G);

			avgError = 0.0;
			value_t N = stabInvG.n_rows * stabInvG.n_rows;
			for (uint_t i = 0; i < stabInvG.n_rows; ++i)
			{
				for (uint_t j = 0; j < stabInvG.n_cols; ++j)
				{
					value_t err = std::abs(invG(i, j) - stabInvG(i, j));
					avgError += err / N;
				}
			}

			invG = stabInvG;
			return 0.0;
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
		
		template<typename Map>
		typename Map::mapped_type& GetWithDef(Map& map, typename Map::key_type key, typename Map::mapped_type defval)
		{
			auto it = map.find( key );
			if (it == map.end())
				map[key] = defval;
			return map[key];
		}
		
	private:
		ConfigSpace_t& configSpace;
		VertexHandler_t vertexHandler;
		matrix_t invG;
		uint_t maxWorms = 2;
		bool print = true;
};
