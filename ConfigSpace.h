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
#include <sstream>
//#define EIGEN_USE_MKL_ALL
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"
#include "measurements.h"
#include "random.h"
#include "parser.h"
#include "types.h"
#include "LookUpTable.h"
#include "UpdateHandler.h"
#include "VertexHandler.h"
#include "GeometryBase.h"

enum UpdateType {AddVertex, RemoveVertex, Add2Vertices, Remove2Vertices, Add5Vertices, Remove5Vertices, Add8Vertices, Remove8Vertices, ZtoW2, W2toZ, ZtoW4, W4toZ, W2toW4, W4toW2, shiftWorm};
enum StateType {Z, W2, W4};

template<typename T>
struct SortSecond
{
	bool operator() (const std::pair<T,T>& lhs, const std::pair<T,T>& rhs)
	{
		return lhs.second < rhs.second;
	}
};

template <typename T>
string ToString ( T arg )
{
	std::ostringstream ss;
	ss << arg;
	return ss.str();
}

template<typename Geometry, typename RNG, typename Value_t, typename Matrix_t>
class ConfigSpace
{
	public:
		typedef Geometry Geometry_t;
		typedef typename Geometry::int_t uint_t;
		typedef typename Geometry::int_t int_t;
		typedef Value_t value_t;
		typedef Matrix_t matrix_t;
		typedef UpdateHandler< ConfigSpace<Geometry, RNG, value_t, Matrix_t> > UpdateHandler_t;
		
		ConfigSpace(RNG& rng)
			:rng(rng), updateHandler(UpdateHandler_t(*this))
		{}

		~ConfigSpace()
		{
			delete lattice;
		}
		
		template<int_t N>
		bool AddRandomVertices(value_t preFactor, bool isWorm)
		{
			if (isWorm)
			{
				updateHandler.GetVertexHandler().template AddRandomWormsToBuffer<N>(nhoodDist);
				preFactor *= updateHandler.GetVertexHandler().VertexBufferParity();
			}
			else
				updateHandler.GetVertexHandler().template AddRandomVerticesToBuffer<N>();
			return updateHandler.template AddVertices<N>(preFactor, isWorm);
		}
		
		template<int_t N>
		bool RemoveRandomVertices(value_t preFactor, bool isWorm)
		{
			if (isWorm)
			{
				if (updateHandler.GetVertexHandler().Worms() < N)
					return false;
				else
				{
					updateHandler.GetVertexHandler().template AddRandomWormIndicesToBuffer<N>();
					preFactor *= updateHandler.GetVertexHandler().WormIndexBufferParity();
					if (updateHandler.GetVertexHandler().template WormIndexBufferDistance<N>())
						return updateHandler.template RemoveVertices<N>(preFactor, isWorm);
					else
						return false;
				}
			}
			else
			{
				if (updateHandler.GetVertexHandler().Vertices() < N)
					return false;
				else
					updateHandler.GetVertexHandler().template AddRandomIndicesToBuffer<N>();
			}
			return updateHandler.template RemoveVertices<N>(preFactor, isWorm);
		}

		template<int_t N>
		bool OpenUpdate()
		{
			if ((state == StateType::Z) && (N == 1) && (updateHandler.GetVertexHandler().Vertices() > 0))
			{
				updateHandler.GetVertexHandler().template AddRandomIndicesToBuffer<N>();
				value_t preFactor = 2.0 * zeta2 * updateHandler.GetVertexHandler().Vertices() / V;
				return updateHandler.OpenUpdate(preFactor);
			}
			else
				return false;
		}

		template<int_t N>
		bool CloseUpdate()
		{
			if ((state == StateType::W2) && (N == 1) && (updateHandler.GetVertexHandler().WormDistance() == 1))
			{
				updateHandler.GetVertexHandler().template AddRandomWormIndicesToBuffer<N>();
				value_t preFactor = V / (2.0 * zeta2 * (updateHandler.GetVertexHandler().Vertices() + 1.0));
				return updateHandler.CloseUpdate(preFactor);
			}
			else
				return false;
		}
		
		template<int_t W>
		bool ShiftWorm()
		{
			return updateHandler.template ShiftWorm<W>();
		}

		/*
		template<int_t W>
		bool ShiftWorm()
		{
			updateHandler.GetVertexHandler().template AddRandomWormIndicesToBuffer<W>();
			updateHandler.GetVertexHandler().ShiftWormToBuffer();
			uint_t m = lattice->MaxDistance();
			value_t preFactorRemove = updateHandler.GetVertexHandler().WormIndexBufferParity() / (lattice->Sites() * m * beta * zeta2);
			value_t preFactorAdd = updateHandler.GetVertexHandler().VertexBufferParity() * lattice->Sites() * m * beta * zeta2;
			value_t detRemove, detAdd;
			updateHandler.template RemoveVertices<W>(preFactorRemove, true, detRemove, UpdateFlag::NoUpdate);
			updateHandler.template AddVertices<W>(preFactorAdd, true, detAdd, UpdateFlag::NoUpdate);
			value_t detShift = std::min({preFactorRemove*detRemove, 1.0}) * std::min({preFactorAdd*detAdd, 1.0});
			//value_t detShift = detRemove*detAdd;

			if (rng() < detShift)
			{
				updateHandler.template RemoveVertices<W>(preFactorRemove, true, detRemove, UpdateFlag::ForceUpdate);
				updateHandler.template AddVertices<W>(preFactorAdd, true, detAdd, UpdateFlag::ForceUpdate);
				return true;
			}
			else
				return false;
		}
		*/

		void Clear()
		{
			updateHandler.Clear();
			state = StateType::Z;
		}
		
		void PrintMatrix(const matrix_t& m)
		{
			for (uint_t i = 0; i < m.rows(); ++i)
			{
				for (uint_t j = 0; j < m.cols(); ++j)
					std::cout << m(i, j) << " ";
				std::cout << std::endl;
			}
		}
		
		StateType State() const
		{
			return state;
		}
		
		/*
		inline value_t LookUpG0(uint_t i1, uint_t i2, value_t tau)
		{
			value_t tau_p;
			if (std::abs(tau) > beta/2.0)
				tau_p = beta - std::abs(tau);
			else
				tau_p = std::abs(tau);
			uint_t t = static_cast<uint_t>(std::abs(tau_p) / dtau);
			uint_t i = std::min({i1, i2}), j = std::max({i1, i2});
			uint_t x = i * lattice->Sites() - (i + i*i) / 2 + j;
			value_t tau_t = t * dtau, G_t = lookUpTableG0(x, t), G_tt = lookUpTableG0(x, t+1);
			value_t g = G_t + (tau_p - tau_t) * (G_tt - G_t) / dtau;
			value_t sign = 1.0;
			bool sameSublattice = lattice->Sublattice(i1) == lattice->Sublattice(i2);
			if (std::abs(tau) > beta/2.0 && (!sameSublattice))
				sign *= -1.0;
			if (tau < 0.0 && sameSublattice)
				sign *= -1.0;
			return sign * g;
		}
		*/
		
		
		inline value_t LookUpG0(uint_t i1, uint_t i2, value_t tau)
		{
			value_t tau_p;
			if (std::abs(tau) > beta/2.0)
				tau_p = beta - std::abs(tau);
			else
				tau_p = std::abs(tau);
			uint_t t = static_cast<uint_t>(std::abs(tau_p) / dtau);
			uint_t x = lookUpIndex(i1, i2);
			value_t tau_t = t * dtau, G_t = lookUpTableG0(x, t), G_tt = lookUpTableG0(x, t+1);
			value_t g = G_t + (tau_p - tau_t) * (G_tt - G_t) / dtau;
			value_t sign = 1.0;
			bool sameSublattice = lattice->Sublattice(i1) == lattice->Sublattice(i2);
			if (std::abs(tau) > beta/2.0 && (!sameSublattice))
				sign *= -1.0;
			if (tau < 0.0 && sameSublattice)
				sign *= -1.0;
			return sign * g;
		}
		

		void EvaluateG0(value_t tau, matrix_t& g0)
		{
			for (uint_t i = 0; i < hopDiag.cols(); ++i)
			{
				value_t ev = evSolver.eigenvalues()[i];
				hopDiag(i, i) = std::exp(-tau * ev) / (1.0 + std::exp(-beta * ev));
			}
			g0 = hopEV * hopDiag * hopEVT;
		}
		
		/*
		void BuildG0LookUpTable()
		{
			uint_t N = lattice->Sites();
			lookUpTableG0.AllocateTable((N*N + N) / 2, nTimeBins + 1);
			matrix_t G0(hopDiag.rows(), hopDiag.cols());
			for (uint_t t = 0; t <= nTimeBins; ++t)
			{
				EvaluateG0(dtau * t, G0);
				for (uint_t i = 0; i < lattice->Sites(); ++i)
				{
					for (uint_t j = i; j < lattice->Sites(); ++j)
					{
						uint_t x = i * N - (i + i*i) / 2 + j;
						lookUpTableG0(x, t) = G0(i, j);
					}
				}
				if (t % (nTimeBins / 10) == 0)
				{
					std::cout << ".";
					std::cout.flush();
				}
			}
		}
		*/
		
		void BuildG0LookUpTable()
		{
			uint_t nValues = BuildG0IndexTable();
			matrix_t G0(hopDiag.rows(), hopDiag.cols());
			value_t threshold = std::pow(10.0, -12.0);
			for (uint_t t = 0; t <= nTimeBins; ++t)
			{
				EvaluateG0(dtau * t, G0);
				std::vector<value_t> cnt(nValues, 0.0);
				for (uint_t i = 0; i < lattice->Sites(); ++i)
				{
					for (uint_t j = i; j < lattice->Sites(); ++j)
					{
						uint_t index = lookUpIndex(i, j);
						lookUpTableG0(index, t) = (G0(i, j) + cnt[index] * lookUpTableG0(index, t)) / (1.0 + cnt[index]);
						//lookUpTableG0(index, t) = G0(i, j);
						++cnt[index];
					}
				}
				if (t % (nTimeBins / 10) == 0)
				{
					std::cout << ".";
					std::cout.flush();
				}
			}
		}
		
		uint_t BuildG0IndexTable()
		{
			uint_t N = lattice->Sites();
			lookUpIndex.AllocateTable(N, N, -1);
			matrix_t G0(hopDiag.rows(), hopDiag.cols());
			value_t threshold = std::pow(10.0, -12.0);
			uint_t t = nTimeBins / 3;
			std::vector<value_t> nValues;
			
			EvaluateG0(dtau * t, G0);
			for (uint_t i = 0; i < lattice->Sites(); ++i)
			{
				for (uint_t j = i; j < lattice->Sites(); ++j)
				{
					for (uint_t k = 0; k < 3; ++k)
					{
						bool isStored = false;
						for (uint_t k = 0; k < nValues.size(); ++k)
						{
							if (std::abs(nValues[k] - G0(i, j)) < threshold)
							{
								isStored = true;
								lookUpIndex(i, j) = k;
								lookUpIndex(j, i) = k;
							}
						}
						if (!isStored)
						{
							nValues.push_back(G0(i, j));
							lookUpIndex(i, j) = nValues.size() - 1;
							lookUpIndex(j, i) = nValues.size() - 1;
						}
					}
				}
			}
			lookUpTableG0.AllocateTable(nValues.size(), nTimeBins + 1, 0.0);
			return nValues.size();
		}

		void BuildHoppingMatrix()
		{
			for (uint_t i = 0; i < lattice->Sites(); ++i)
			{
				for (uint_t j = 0; j < lattice->Sites(); ++j)
				{
					if (lattice->IsNeighbor(i, j))
						hoppingMatrix(i, j) = -t;
					else
						hoppingMatrix(i, j) = 0;
				}
			}
			evSolver.compute(hoppingMatrix);
			hopDiag = matrix_t::Zero(hoppingMatrix.rows(), hoppingMatrix.cols());
			hopEV = evSolver.eigenvectors();
			hopEVT = evSolver.eigenvectors().adjoint();
		}

		void ResizeGeometry(uint_t l)
		{
			L = l;
			lattice->Resize(l, rng);
			hoppingMatrix.resize(lattice->Sites(), lattice->Sites());
		}
		
		void SetTemperature(value_t T)
		{
			beta = 1.0 / T;
			dtau = beta / (2 * static_cast<value_t>(nTimeBins));
			infinTau = dtau / 1000.0;
		}
		
		void Serialize(odump& d)
		{
			d.write(state);
			updateHandler.Serialize(d);
		}
		
		void Serialize(idump& d)
		{
			d.read(state);
			updateHandler.Serialize(d);
		}
		
		value_t AdditionFactorialRatio(uint_t k, uint_t n)
		{
			if (k <= 0 || n <= 0)
				return 1.0;
			value_t result = 1.0;
			for (uint_t i = 1; i <= n; ++i)
				result /= static_cast<value_t>(k + i);
			return result;
		}

		value_t RemovalFactorialRatio(uint_t k, uint_t n)
		{
			if (k <= 0 || n <= 0)
				return 1.0;
			value_t result = 1.0;
			for (uint_t i = 1; i <= n; ++i)
				result *= k - n + i;
			return result;
		}
	public:
		RNG& rng;
		uint_t L;
		Geometry* lattice;
		UpdateHandler_t updateHandler;
		value_t beta;
		value_t V;
		value_t t;
		value_t infinTau;
		value_t zeta2;
		value_t zeta4;
		uint_t nTimeBins;
		uint_t maxWorms = 4;
		//LookUpTable<value_t, 2> lookUpTableG0;
		LookUpTable<value_t, 2> lookUpTableG0;
		LookUpTable<int_t, 2> lookUpIndex;
		value_t dtau;
		StateType state = StateType::Z;
		matrix_t hoppingMatrix;
		matrix_t hopDiag;
		matrix_t hopEV;
		matrix_t hopEVT;
		Eigen::SelfAdjointEigenSolver<matrix_t> evSolver;
		uint_t nhoodDist;
		bool fileIO;
};
