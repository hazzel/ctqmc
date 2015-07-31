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
		{
			invG.resize(200, 200);
			G.resize(200, 200);
		}
		
		void Init()
		{}
		
		template<typename Matrix_t>
		void PrintMatrix(const Matrix_t& M)
		{
			for (uint_t i = 0; i < M.n_rows; ++i)
			{
				for (uint_t j = 0; j < M.n_cols; ++j)
					std::cout << M(i, j) << " ";
				std::cout << std::endl;
			}
		}

		template<int_t N>
		bool AddVertices(value_t preFactor, bool isWorm, bool force = false)
		{
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;
			
			if (k + n > invG.n_rows)
			{
				matrix_t buf = invG;
				invG.resize(k + n + 10, k + n + 10);
				invG.submat(0, 0, buf.n_rows - 1, buf.n_cols - 1) = buf;
				buf = G;
				G.resize(k + n + 10, k + n + 10);
				G.submat(0, 0, buf.n_rows - 1, buf.n_cols - 1) = buf;
			}
			matrix_t u(k, n), v(n, k), a(n, n);
			vertexHandler.WoodburyAddVertices(u, v, a);

			matrix_t invGu(k, n);
			matrix_t invS(n, n);
			if (k > 0)
			{
				invGu = invG.submat(0, 0, k - 1, k - 1) * u;
				invS = a - v * invGu;
			}
			else
				invS = a;

			value_t acceptRatio = preFactor * arma::det(invS);
			if (print && acceptRatio < 0.0)
			{
				std::cout << "AddVertices(" << N << "): AcceptRatio: " << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
				StabilizeInvG();
			}
			if (force)
				acceptRatio = 1.0;
			if (configSpace.rng() < acceptRatio)
			{
				matrix_t S = arma::inv(invS);
				if (k > 0)
				{
					matrix_t vinvG = invG.submat(0, 0, k - 1, k - 1).t() * v.t();
					arma::inplace_trans(vinvG);
					invG.submat(k, 0, n + k - 1, k - 1) = -S * vinvG;
					invG.submat(0, 0, k - 1, k - 1) = invG.submat(0, 0, k - 1, k - 1) - invGu * invG.submat(k, 0, k + n - 1, k - 1);
					invG.submat(0, k, k - 1, k + n - 1) = -invGu * S;
				}
				invG.submat(k, k, k + n - 1, k + n - 1) = S;
				if (k > 0)
				{
					G.submat(0, k, k - 1, k + n - 1) = u;
					G.submat(k, 0, k + n - 1, k - 1) = v;
				}
				G.submat(k, k, k + n - 1, k + n - 1) = a;
				
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

			matrix_t S(n, n);
			vertexHandler.FillSMatrix(S, invG, isWorm);
			value_t acceptRatio = preFactor * arma::det(S);
			if (print && acceptRatio < 0.0)
			{
				std::cout << "RemoveVertices(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
				StabilizeInvG();
			}
			if (configSpace.rng() < acceptRatio)
			{
				if (k != n)
				{
					vertexHandler.PermuteProgagatorMatrix(invG, isWorm);
					
					arma::inplace_trans(S);
					matrix_t t = invG.submat(k - n, 0, k - 1, k - n - 1).t() * arma::inv(S);
					arma::inplace_trans(t);
					invG.submat(0, 0, k - n - 1, k - n - 1) -= invG.submat(0, k - n, k - n - 1, k - 1) * t;
					
					vertexHandler.PermuteProgagatorMatrix(G, isWorm);
				}
				vertexHandler.RemoveBufferedVertices2(isWorm);
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

			matrix_t shiftedWormU(k, l);
			matrix_t shiftedWormV(l, k);
			matrix_t shiftedWormA(l, l);

			vertexHandler.ShiftWormToBuffer();
			vertexHandler.PermuteProgagatorMatrix(invG, true);
			vertexHandler.PermuteProgagatorMatrix(G, true);
			vertexHandler.PermuteVertices(true);
			vertexHandler.WoodburyShiftWorm(shiftedWormU, shiftedWormV, shiftedWormA);

			matrix_t t = invG.submat(k, 0, k + l - 1, k - 1).t() * arma::inv(invG.submat(k, k, k + l - 1, k + l - 1)).t();
			arma::inplace_trans(t);
			matrix_t M = invG.submat(0, 0, k - 1, k - 1) - invG.submat(0, k, k - 1, k + l - 1) * t;
			matrix_t Mt = M.t();
			
			arma::inplace_trans(shiftedWormV);
			matrix_t shiftedInvS = shiftedWormA - (shiftedWormU.t() * Mt * shiftedWormV).t();
			value_t detShiftedInvS = arma::det(shiftedInvS);

			value_t acceptRatio = detShiftedInvS * arma::det(invG.submat(k, k, k + l - 1, k + l - 1)) * vertexHandler.WormShiftParity();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "WormShift: AcceptRatio: " << acceptRatio << std::endl;
				StabilizeInvG();
			}
			//if (configSpace.rng() < std::abs(p1) && configSpace.rng() < std::abs(p2))
			if (configSpace.rng() < acceptRatio)
			{
				matrix_t S = arma::inv(shiftedInvS);
				matrix_t vM = Mt * shiftedWormV;
				arma::inplace_trans(vM);
				matrix_t Mu = M * shiftedWormU;
				invG.submat(k, 0, k + l - 1, k - 1) = -S * vM;
				invG.submat(0, 0, k - 1, k - 1) = M - Mu * invG.submat(k, 0, k + l - 1, k - 1);
				invG.submat(0, k, k - 1, k + l - 1) = -Mu * S;
				invG.submat(k, k, k + l - 1, k + l - 1) = S;
				
				G.submat(0, k, k - 1, k + l - 1) = shiftedWormU;
				G.submat(k, 0, k + l - 1, k - 1) = shiftedWormV.t();
				G.submat(k, k, k + l - 1, k + l - 1) = shiftedWormA;
				
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
			//preFactorRem = 1.0;
			//preFactorAdd = 1.0;
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
		
		value_t MeasureM2()
		{
			return 0.;
		}
		
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
			matrix_t t = G.submat(0, 0, k - 1, k - 1) * invG.submat(0, 0, k - 1, k - 1);
			uint_t i = static_cast<uint_t>(configSpace.rng() * (t.n_rows - 1)), j = static_cast<uint_t>(configSpace.rng()  * (t.n_cols - 1));
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
			invG.submat(0, 0, k - 1, k - 1) = arma::inv(G.submat(0, 0, k - 1, k - 1));
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
			matrix_t g(k, k);
			vertexHandler.PropagatorMatrix(g);
			if (k > G.n_rows)
			{
				G.resize(k, k);
				invG.resize(k, k);
			}
			if (k > 0)
			{
				G.submat(0, 0, k - 1, k - 1) = g;
				StabilizeInvG();
			}
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
		matrix_t G;
		uint_t matSize = 0;
		uint_t maxWorms = 2;
		bool print = true;
};
