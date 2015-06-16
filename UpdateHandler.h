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
			for (uint_t i = 0; i < M.n_rows; ++i)
			{
				for (uint_t j = 0; j < M.n_cols; ++j)
					std::cout << M(i, j) << " ";
				std::cout << std::endl;
			}
		}

		template<int_t N>
		bool AddVertices(value_t preFactor, bool isWorm)
		{
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;
			matrix_t u(k, n), v(n, k), a(n, n);
			vertexHandler.WoodburyAddVertices(u, v, a);

			matrix_t invGu = invG * u;
			matrix_t invS = a - v * invGu;

			value_t acceptRatio = preFactor * arma::det(invS);
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
			if (isWorm && vertexHandler.Worms() < N)
				return false;
			if ((!isWorm) && vertexHandler.Vertices() < N)
				return false;
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;

			matrix_t S(n, n);
			vertexHandler.FillSMatrix(S, invG, isWorm);
			value_t acceptRatio = preFactor * arma::det(S);
			if (print && acceptRatio < 0.0)
			{
				std::cout << "RemoveVertices(" << N << "): AcceptRatio" << acceptRatio << std::endl;
				std::cout << "IsWorm: " << isWorm << ", Vertices: " << vertexHandler.Vertices() << ", Worms: " << vertexHandler.Worms() << std::endl;
			}
			if (configSpace.rng() < acceptRatio)
			{
				arma::vec perm(k);
				vertexHandler.PermutationMatrix(perm, isWorm);
				matrix_t invGp(k, k);
				
				for (uint_t i = 0; i < k; ++i)
					for (uint_t j = 0; j < k; ++j)
						invGp(i, j) = invG(perm(i), perm(j));
				/*
				std::vector<std::size_t> pos;
				if (vertexHandler.IndexBuffer()[0] > 0)
					pos.push_back(0);
				for (uint_t i = 0; i < N; ++i)
				{
					if (vertexHandler.IndexBuffer()[0] > 0)
						pos.push_back(vertexHandler.IndexBuffer()[2*i] - 1);
					if (vertexHandler.IndexBuffer()[2*N - 2] < k - 2)
						pos.push_back(vertexHandler.IndexBuffer()[2*i] + 2);
				}
				if (vertexHandler.IndexBuffer()[2*N - 2] < k - 2)
					pos.push_back(k - 1);
				for (uint_t i = 0; i < pos.size(); i+=2)
				{
					for (uint_t j = 0; j < pos.size(); j+=2)
					{
						invGp.submat(pos[i], pos[j], pos[i+1], pos[j+1]);
					}
				}
				*/
				
				matrix_t t = arma::inv(S) * invGp.submat(k - n, 0, k - 1, k - n - 1);
				invG = invGp.submat(0, 0, k - n - 1, k - n - 1) - invGp.submat(0, k - n, k - n - 1, k - 1) * t;

				vertexHandler.RemoveBufferedVertices(isWorm);
				if (k - n > 4)
				{
					matrix_t cols = invG.submat(0, 0, k - n - 1, 1);
					invG.submat(0, 0, k - n - 1, 1) = invG.submat(0, 2, k - n - 1, 3);
					invG.submat(0, 2, k - n - 1, 3) = cols;
					matrix_t rows = invG.submat(0, 0, 1, k - n - 1);
					invG.submat(0, 0, 1, k - n - 1) = invG.submat(2, 0, 3, k - n - 1);
					invG.submat(2, 0, 3, k - n - 1) = rows;
					vertexHandler.SwapVertexPosition(0, 2);
				}
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
			const uint_t l = 2 * W;
			
			vertexHandler.ShiftWormToBuffer();
			matrix_t U_1(k + l, l);
			matrix_t V_2(l, k + l);
			vertexHandler.WoodburyShiftRowsCols(U_1, V_2);
			matrix_t I = arma::eye<matrix_t>(l, l);
			matrix_t invGu1 = invG * U_1;
			matrix_t invGv1(l, k + l);
			matrix_t invGu2(k + l, l);
			for (uint_t i = 0; i < l; i+=2)
			{
				uint_t pos = vertexHandler.WormPositions()[i];
				invGv1.submat(i, 0, i + 1, k + l - 1) = invG.submat(pos, 0, pos + 1, k + l - 1);
				invGu2.submat(0, i, k + l - 1, i + 1) = invG.submat(0, pos, k + l - 1, pos + 1);
			}
			
			matrix_t w1 = I + invGv1 * U_1;
			matrix_t w2 = I + V_2 * invGu2;
			matrix_t w3(l, l);
			for (uint_t i = 0; i < l; i+=2)
			{
				uint_t pos = vertexHandler.WormPositions()[i];
				w3.submat(i, 0, i + 1, l - 1) = invGu2.submat(pos, 0, pos + 1, l - 1);
			}
			matrix_t w4 = invGu1 * arma::inv(w1);
			matrix_t w5 = V_2 * w4 * w3;
			
			value_t acceptRatio = arma::det(w1) * arma::det(w2 - w5) * vertexHandler.WormShiftParity();
			if (print && acceptRatio < 0.0)
			{
				std::cout << "WormShift: AcceptRatio: " << acceptRatio << std::endl;
			}
			if (configSpace.rng() < acceptRatio)
			{
				//update columns
				invG -= w4 * invGv1;
				
				//update rows
				matrix_t invGv2 = V_2 * invG;
				matrix_t w6(l, l);
				for (uint_t i = 0; i < l; i+=2)
				{
					uint_t pos = vertexHandler.WormPositions()[i];
					invGu2.submat(0, i, k + l - 1, i + 1) = invG.submat(0, pos, k + l - 1, pos + 1);
					w6.submat(0, i, l - 1, i + 1) = invGv2.submat(0, pos, l - 1, pos + 1);
				}
				invG -= invGu2 * arma::inv(I + w6) * invGv2;
				vertexHandler.ApplyWormShift();
				
				return true;
			}
			else
			{
				return false;
			}
		}
		
		template<int_t W>
		bool ReplaceWorm()
		{
			return vertexHandler.ReplaceWorm<W>();
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
					value_t det;
					AddVertices<W>(preFactorAdd * vertexHandler.VertexBufferParity(), true, det, UpdateFlag::ForceUpdate);
					return false;
				}
			}
			return false;
		}
		*/
		
		value_t GetWeight()
		{
			return 1.0;
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
