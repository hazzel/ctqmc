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
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"
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
		template<int_t N, int_t M> using matrix_t = Eigen::Matrix<value_t, N, M, Eigen::RowMajor>;
		template<int_t N> using inv_solver_t = Eigen::FullPivLU< matrix_t<N, N> >;
		
		UpdateHandler(ConfigSpace_t& configSpace)
			: configSpace(configSpace), vertexHandler(VertexHandler_t(configSpace))
		{
			Ubuf.resize(100, 100);
			Vbuf.resize(100, 100);
			Abuf.resize(20, 20);
		}
		
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
		
		/*
		template<int_t N>
		bool AddVertices(value_t preFactor, bool isWorm, value_t& det, UpdateFlag flag)
		{
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;
			if (k > Ubuf.rows() || k > Vbuf.cols())
			{
				Ubuf.resize(k+20, Ubuf.cols());
				Vbuf.resize(Vbuf.rows(), k+20);
			}
			matrix_t<n, n> a(n, n);
			vertexHandler.WoodburyAddVertices(Ubuf, Vbuf, a);

			matrix_t<Eigen::Dynamic, n> invGu = invG * Ubuf.topLeftCorner(k, n);
			matrix_t<n, n> invS = a.template topLeftCorner<n, n>();
			invS.noalias() -= Vbuf.topLeftCorner(n, k) * invGu;

			value_t acceptRatio;
			if (flag == UpdateFlag::NormalUpdate)
			{
				det = invS.determinant();
				acceptRatio = preFactor * det;
			}
			else if (flag == UpdateFlag::NoUpdate)
			{
				det = invS.determinant();
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
				matrix_t<n, n> S = invS.inverse();
				matrix_t<n, Eigen::Dynamic> vinvG = Vbuf.topLeftCorner(n, k) * invG;
				
				invG.conservativeResize(k + n, k + n);
				invG.bottomLeftCorner(n, k) = -S * vinvG;
				invG.topLeftCorner(k, k).noalias() -= invGu * invG.bottomLeftCorner(n, k);
				invG.topRightCorner(k, n).noalias() = -invGu * S;
				invG.template bottomRightCorner<n, n>() = S;
				
				vertexHandler.AddBufferedVertices(isWorm);
				return true;
			}
			return false;
		}
		*/
		
		template<int_t N>
		bool AddVertices(value_t preFactor, bool isWorm, value_t& det, UpdateFlag flag)
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

			value_t acceptRatio;
			if (flag == UpdateFlag::NormalUpdate)
			{
				det = invS.determinant();
				acceptRatio = preFactor * det;
			}
			else if (flag == UpdateFlag::NoUpdate)
			{
				det = invS.determinant();
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
				matrix_t<n, n> S = invS.inverse();
				matrix_t<n, Eigen::Dynamic> vinvG = v * invG;
				
				invG.conservativeResize(k + n, k + n);
				invG.bottomLeftCorner(n, k) = -S * vinvG;
				invG.topLeftCorner(k, k).noalias() -= invGu * invG.bottomLeftCorner(n, k);
				invG.topRightCorner(k, n).noalias() = -invGu * S;
				invG.template bottomRightCorner<n, n>() = S;
				
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
			
			std::vector<value_t> perm(k);
			vertexHandler.PermutationMatrix(perm, isWorm);
			
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> invGp(k, k);
			for (uint_t i = 0; i < k; ++i)
				for (uint_t j = 0; j < k; ++j)
					invGp(i, j) = invG(perm[i], perm[j]);
			
			matrix_t<n, n> S = invGp.template bottomRightCorner<n, n>();
			value_t acceptRatio;
			if (flag == UpdateFlag::NormalUpdate)
			{
				det = S.determinant();
				acceptRatio = preFactor * det;
			}
			else if (flag == UpdateFlag::NoUpdate)
			{
				det = S.determinant();
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
				matrix_t<n, Eigen::Dynamic> t = S.inverse() * invGp.bottomLeftCorner(n, k - n);
				invG = invGp.topLeftCorner(k - n, k - n);
				invG.noalias() -= invGp.topRightCorner(k - n, n) * t;

				vertexHandler.RemoveBufferedVertices(isWorm);
				return true;
			}
			else
			{
				return false;
			}
		}
		
		/*
		template<int_t N>
		bool RemoveVertices(value_t preFactor, bool isWorm, value_t& det, UpdateFlag flag)
		{
			if (isWorm && vertexHandler.Worms() < N)
				return false;
			if ((!isWorm) && vertexHandler.Vertices() < N)
				return false;
			uint_t k = 2 * (vertexHandler.Vertices() + vertexHandler.Worms());
			const uint_t n = 2 * N;
			
			std::vector<value_t> perm(k);
			vertexHandler.PermutationMatrix(perm, isWorm);
			
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> invGp(k, k);
			for (uint_t i = 0; i < k; ++i)
				for (uint_t j = 0; j < k; ++j)
					invGp(i, j) = invG(perm[i], perm[j]);
			
			matrix_t<n, n> S = invGp.template bottomRightCorner<n, n>();
			value_t acceptRatio;
			if (flag == UpdateFlag::NormalUpdate)
			{
				det = S.determinant();
				acceptRatio = preFactor * det;
			}
			else if (flag == UpdateFlag::NoUpdate)
			{
				det = S.determinant();
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
				matrix_t<n, Eigen::Dynamic> t = S.inverse() * invGp.bottomLeftCorner(n, k - n);
				invG = invGp.topLeftCorner(k - n, k - n);
				invG.noalias() -= invGp.topRightCorner(k - n, n) * t;

				vertexHandler.RemoveBufferedVertices(isWorm);
				return true;
			}
			else
			{
				return false;
			}
		}
		*/

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

			std::vector<value_t> perm(k + l);
			vertexHandler.PermutationMatrix(perm, true);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> invGp(k + l, k + l);
			for (uint_t i = 0; i < k + l; ++i)
				for (uint_t j = 0; j < k + l; ++j)
					invGp(i, j) = invG(perm[i], perm[j]);
			
			
			//Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(k + l);
			//vertexHandler.PermutationMatrix(perm.indices(), true);
			//invG = perm.transpose() * invG * perm;
			
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> M = invGp.topLeftCorner(k, k);
			M.noalias() -= invGp.topRightCorner(k, l) * invGp.template bottomRightCorner<l, l>().inverse() * invGp.bottomLeftCorner(l, k);
			
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> g(k + l, k + l);
			vertexHandler.PropagatorMatrix(g);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> gp(k + l, k + l);
			for (uint_t i = 0; i < k + l; ++i)
				for (uint_t j = 0; j < k + l; ++j)
					gp(i, j) = g(perm[i], perm[j]);
			
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> gtl = gp.topLeftCorner(k, k);
			int mc = MatrixCondition(gtl);
			std::cout << mc << std::endl;
			if (mc > 100000)
				std::cin.get();

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
				invGp.topLeftCorner(k, k) = M;
				invGp.topLeftCorner(k, k).noalias() -= Mu * R;
				invGp.topRightCorner(k, l).noalias() = -Mu * S;
				invGp.bottomLeftCorner(l, k) = R;
				invGp.template bottomRightCorner<l, l>() = S;
				//invG = perm * invG * perm.transpose();
				
				for (uint_t i = 0; i < k + l; ++i)
					for (uint_t j = 0; j < k + l; ++j)
						invG(perm[i], perm[j]) = invGp(i, j);

				vertexHandler.ApplyWormShift();
				
				
				//value_t det;
				//RemoveVertices<W>(preFactorRem * vertexHandler.WormIndexBufferParity(), true, det, UpdateFlag::ForceUpdate);
				//AddVertices<W>(preFactorAdd * vertexHandler.VertexBufferParity(), true, det, UpdateFlag::ForceUpdate);
				
				return true;
			}
			else
			{
				//invG = perm * invG * perm.transpose();
				
				for (uint_t i = 0; i < k + l; ++i)
					for (uint_t j = 0; j < k + l; ++j)
						invG(perm[i], perm[j]) = invGp(i, j);
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
			if (M.rows() == 0)
				return 0.0;
			Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
			return svd.singularValues()(0) / svd.singularValues()(M.rows()-1);
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
		
		value_t StabilizeInvG(value_t& avgError)
		{
			if (invG.rows() == 0)
				return 0.0;
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> G(invG.rows(), invG.cols());
			vertexHandler.PropagatorMatrix(G);
			inv_solver_t<Eigen::Dynamic> solver(G);
			matrix_t<Eigen::Dynamic, Eigen::Dynamic> stabInvG = solver.inverse();

			avgError = 0.0;
			value_t N = stabInvG.rows() * stabInvG.rows();
			for (uint_t i = 0; i < stabInvG.rows(); ++i)
			{
				for (uint_t j = 0; j < stabInvG.cols(); ++j)
				{
					value_t err = std::abs(invG(i, j) - stabInvG(i, j));
					avgError += err / N;
				}
			}

			invG = stabInvG;
			//return MatrixCondition(invG);
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
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> invG;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> Ubuf;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> Vbuf;
		matrix_t<Eigen::Dynamic, Eigen::Dynamic> Abuf;
		uint_t maxWorms = 2;
		bool print = true;
};
