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
#include "Eigen/Dense"
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
		template<int_t N, int_t M> using matrix_t = Eigen::Matrix<value_t, N, M, Eigen::ColMajor>;
		template<int_t N> using inv_solver_t = Eigen::FullPivLU< matrix_t<N, N> >;
		using dmatrix_t = matrix_t<Eigen::Dynamic, Eigen::Dynamic>;
		
		UpdateHandler(ConfigSpace_t& configSpace)
			: configSpace(configSpace), vertexHandler(VertexHandler_t(configSpace))
		{
			invG.resize(200, 200);
			G.resize(200, 200);
		}
		
		void Init()
		{}
		
		template<typename Matrix_t>
		void PrintMatrix(const Matrix_t& M)
		{
			for (uint_t i = 0; i < M.rows(); ++i)
			{
				for (uint_t j = 0; j < M.cols(); ++j)
					std::cout << M(i, j) << " ";
				std::cout << std::endl;
			}
		}

		template<int_t N>
		bool AddVertices(value_t preFactor, bool isWorm, bool force = false)
		{
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;
			
			if (k + n > invG.rows())
			{
				dmatrix_t buf = invG;
				invG.resize(k + n + 10, k + n + 10);
				invG.topLeftCorner(buf.rows(), buf.cols()) = buf;
				buf = G;
				G.resize(k + n + 10, k + n + 10);
				G.topLeftCorner(buf.rows(), buf.cols()) = buf;
			}
			dmatrix_t u(k, n), v(n, k);
			matrix_t<n, n> a(n, n);
			vertexHandler.WoodburyAddVertices(u, v, a);

			dmatrix_t invGu(k, n);
			matrix_t<n, n> invS(n, n);
			if (k > 0)
			{
				invGu = invG.topLeftCorner(k, k) * u;
				invS = a - v * invGu;
			}
			else
				invS = a;

			value_t acceptRatio = preFactor * invS.determinant();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "AddVertices(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
			}
			if (force)
				acceptRatio = 1.0;
			if (configSpace.rng() < acceptRatio)
			{
				dmatrix_t S = invS.inverse();
				if (k > 0)
				{
					dmatrix_t vinvG = v * invG.topLeftCorner(k, k);
					invG.block(k, 0, n, k) = -S * vinvG;
					invG.topLeftCorner(k, k) -= invGu * invG.block(k, 0, n, k);
					invG.block(0, k, k, n) = -invGu * S;
				}
				invG.template block<n, n>(k, k) = S;
				if (k > 0)
				{
					G.block(0, k, k, n) = u;
					G.block(k, 0, n, k) = v;
				}
				G.template block<n, n>(k, k) = a;
				
				vertexHandler.AddBufferedVertices(isWorm);
				return true;
			}
			return false;
		}

		template<int_t N>
		bool RemoveVertices(value_t preFactor, bool isWorm)
		{
			if (isWorm && vertexHandler.Worms() < N)
				return false;
			if ((!isWorm) && vertexHandler.Vertices() < N)
				return false;
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;
			if (k == 0)
				return false;

			matrix_t<n, n> S(n, n);
			vertexHandler.FillSMatrix(S, invG, isWorm);
			value_t acceptRatio = preFactor * S.determinant();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "RemoveVertices(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
			}
			if (configSpace.rng() < acceptRatio)
			{
				if (k != n)
				{
					vertexHandler.PermuteProgagatorMatrix(invG, isWorm);
					
					dmatrix_t t = S.inverse() * invG.block(k - n, 0, n, k - n);
					invG.topLeftCorner(k - n, k - n) -= invG.block(0, k - n, k - n, n) * t;
					
					vertexHandler.PermuteProgagatorMatrix(G, isWorm);
				}
				vertexHandler.RemoveBufferedVertices2(isWorm);
				StabilizeInvG();
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
		
		
		template<int_t W>
		bool ShiftWorm()
		{
			uint_t k = 2 * vertexHandler.Vertices();
			if (k == 0)
				return false;
			const uint_t l = 2 * W;

			dmatrix_t shiftedWormU(k, l);
			dmatrix_t shiftedWormV(l, k);
			dmatrix_t shiftedWormA(l, l);

			vertexHandler.ShiftWormToBuffer();
			vertexHandler.PermuteProgagatorMatrix(invG, true);
			vertexHandler.PermuteProgagatorMatrix(G, true);
			vertexHandler.PermuteVertices(true);
			vertexHandler.WoodburyShiftWorm(shiftedWormU, shiftedWormV, shiftedWormA);

			dmatrix_t t = invG.template block<l, l>(k, k).inverse() * invG.block(k, 0, l, k);
			dmatrix_t M = invG.topLeftCorner(k, k) - invG.block(0, k, k, l) * t;
			
			dmatrix_t shiftedInvS = shiftedWormA - shiftedWormV * M * shiftedWormU;
			value_t detShiftedInvS = shiftedInvS.determinant();

			value_t acceptRatio = detShiftedInvS * invG.template block<l, l>(k, k).determinant() * vertexHandler.WormShiftParity();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "WormShift: AcceptRatio: " << acceptRatio << std::endl;
			}
			if (configSpace.rng() < acceptRatio)
			{
				matrix_t<l, l> S = shiftedInvS.inverse();
				dmatrix_t vM = shiftedWormV * M;
				dmatrix_t Mu = M * shiftedWormU;
				invG.block(k, 0, l, k) = -S * vM;
				invG.topLeftCorner(k, k) = M - Mu * invG.block(k, 0, l, k);
				invG.block(0, k, k, l) = -Mu * S;
				invG.template block<l, l>(k, k) = S;
				
				G.block(0, k, k, l) = shiftedWormU;
				G.block(k, 0, l, k) = shiftedWormV;
				G.template block<l, l>(k, k) = shiftedWormA;
				
				vertexHandler.ApplyWormShift();
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
			preFactorRem = 1.0;
			preFactorAdd = 1.0;
			if (RemoveVertices<W>(preFactorRem * vertexHandler.WormIndexBufferParity(), true))
			{
				if (AddVertices<W>(preFactorAdd * vertexHandler.VertexBufferParity(), true))
				{
					return true;
				}
				else
				{
					vertexHandler.template RestoreAfterShift<W>();
					AddVertices<W>(preFactorAdd * vertexHandler.VertexBufferParity(), true, true);
					return false;
				}
			}
			return false;
		}
		*/
		
		template<int_t W>
		bool ReplaceWorm()
		{
			return vertexHandler.ReplaceWorm<W>();
		}
		
		value_t GetWeight()
		{
			return 1.0;
		}
		
		void Clear()
		{
			vertexHandler.Clear();
		}
		
		template<typename Matrix>
		value_t MatrixCondition(Matrix& M)
		{
			return 0.0;
		}
		
		value_t IsStableInverse()
		{
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			if (k == 0)
				return 0.0;
			dmatrix_t t = G.topLeftCorner(k, k) * invG.topLeftCorner(k, k);
			uint_t i = static_cast<uint_t>(configSpace.rng() * t.rows()), j = static_cast<uint_t>(configSpace.rng()  * t.cols());
			value_t ref;
			if (i == j)
				ref = 1.0;
			else
				ref = 0.0;
			value_t err = std::abs(t(i, j) - ref);
			if (err > std::pow(10.0, -8.0))
			{
				std::cout << "Warning! Round off error at (" << i << ", " << j << ") of " << err << std::endl;
				std::cout << "Stabilization necessary." << std::endl;
				StabilizeInvG();
			}
			else if (err > std::pow(10.0, -10.0))
			{
				std::cout << "Warning! Stabilization necessary." << std::endl;
				StabilizeInvG();
			}
			return err;
		}

		value_t StabilizeInvG()
		{
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			if (k == 0)
				return 0.0;
			invG.topLeftCorner(k, k) = G.topLeftCorner(k, k).inverse();
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
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			dmatrix_t g(k, k);
			vertexHandler.PropagatorMatrix(g);
			if (k > G.rows())
			{
				G.resize(k, k);
				invG.resize(k, k);
			}
			G.topLeftCorner(k, k) = g;
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
		dmatrix_t invG;
		dmatrix_t G;
		uint_t matSize = 0;
		uint_t maxWorms = 2;
		bool print = true;
};
