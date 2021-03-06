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

//#include <unsupported/Eigen/MPRealSupport>
//using namespace mpfr;

enum UpdateType {AddVertex, RemoveVertex, AddTwoVertices, RemoveTwoVertices, ZtoW2, W2toZ, ZtoW4, W4toZ, W2toW4, W4toW2, shiftWorm};
enum StateType {Z, W2, W4};

template<typename T>
struct SortSecond
{
	bool operator() (const std::pair<T,T>& lhs, const std::pair<T,T>& rhs)
	{
		return lhs.second < rhs.second;
	}
};

template<typename Geometry, typename RNG, typename Value_t, typename Matrix_t>
class ConfigSpace
{
	public:
		typedef Geometry Geometry_t;
		typedef typename Geometry::index_t uint_t;
		typedef typename Geometry::int_t int_t;
		typedef Value_t value_t;
		typedef Matrix_t matrix_t;
		typedef UpdateHandler< ConfigSpace<Geometry, RNG, value_t, Matrix_t> > UpdateHandler_t;
		
		ConfigSpace(RNG& rng)
			:rng(rng), updateHandler(UpdateHandler_t(*this))
		{
			updateList.resize(10, "");
		}
		
		template<int_t N>
		bool AddRandomVertices()
		{
			updateHandler.GetVertexHandler().template AddRandomVerticesToBuffer<N>();
			if (updateHandler.GetVertexHandler().Worms() == 0)
			{
				updateList.pop_front();
				updateList.push_back("AddVertices N=" + std::to_string(N));
				return updateHandler.template AddVertices<N>();
			}
			else
			{
				updateList.pop_front();
				updateList.push_back("AddVerticesWithWorms N=" + std::to_string(N));
				return updateHandler.template AddVerticesWithWorms<N>();
			}
		}
		
		template<int_t N>
		bool RemoveRandomVertices()
		{
			if (updateHandler.GetVertexHandler().Vertices() < N)
				return false;
			updateHandler.GetVertexHandler().template AddRandomIndicesToBuffer<N>();
			if (updateHandler.GetVertexHandler().Worms() == 0)
			{
				updateList.pop_front();
				updateList.push_back("RemoveVertices N=" + std::to_string(N));
				return updateHandler.template RemoveVertices<N>();
			}
			else
			{
				updateList.pop_front();
				updateList.push_back("RemoveVerticesWithWorms N=" + std::to_string(N));
				return updateHandler.template RemoveVerticesWithWorms<N>();
			}
		}
		
		template<int_t N>
		bool AddRandomWorms(value_t preFactor)
		{
			updateList.pop_front();
			updateList.push_back("AddWorms N=" + std::to_string(N));
			updateHandler.GetVertexHandler().template AddRandomWormsToBuffer<N>();
			return updateHandler.template AddWorms<N>(preFactor);
		}
		
		template<int_t N>
		bool RemoveRandomWorms(value_t preFactor)
		{
			if (updateHandler.GetVertexHandler().Worms() < N)
				return false;
			updateList.pop_front();
			updateList.push_back("RemoveWorms N=" + std::to_string(N));
			updateHandler.GetVertexHandler().template AddRandomWormIndicesToBuffer<N>();
			return updateHandler.template RemoveWorms<N>(preFactor);
		}
		
		bool ShiftWorm()
		{
			updateList.pop_front();
			updateList.push_back("ShiftWorm");
			return updateHandler.ShiftWorm();
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
		
		void PrintLastUpdates()
		{
			std::cout << "Last Updates:" << std::endl;
			for (auto str : updateList)
				std::cout << str << std::endl;
			std::cout << std::endl;
		}

		StateType State() const
		{
			return state;
		}
		
		value_t LookUpG0(uint_t i1, uint_t i2, value_t tau)
		{
			uint_t i = static_cast<uint_t>(std::abs(tau) / dtau);
			value_t tau_i = i * dtau;
			uint_t dist = lattice.Distance(i1, i2);
			value_t g = lookUpTableG0[dist][i] + (std::abs(tau) - tau_i) * lookUpTableDtG0[dist][i];
			if (tau >= 0.0)
			{
				return g;
			}
			else
			{
				if (lattice.Sublattice(i1) == lattice.Sublattice(i2))
					return -g;
				else
					return g;
			}
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
		
		void BuildG0LookUpTable()
		{
			uint_t i = lattice.RandomSite(rng);
			std::vector<uint_t> sites;
			for (uint_t r = 0; r <= lattice.MaxDistance(); ++r)
			{
				for (int_t j = 0; j < lattice.Sites(); ++j)
				{
					if (lattice.Distance(i, j) == r)
					{
						sites.push_back(j);
						break;
					}
				}
			}
			matrix_t G0(hopDiag.rows(), hopDiag.cols());
			for (uint_t t = 0; t <= nTimeBins; ++t)
			{
				EvaluateG0(dtau * t, G0);
				for (uint_t r = 0; r < sites.size(); ++r)
					lookUpTableG0[r][t] = G0(i, sites[r]);
				if (t % (nTimeBins / 3) == 0)
				{
					std::cout << ".";
					std::cout.flush();
				}
			}

			for (uint_t t = 0; t < nTimeBins; ++t)
				for (uint_t r = 0; r < sites.size(); ++r)
					lookUpTableDtG0[r][t] = (lookUpTableG0[r][t + 1] - lookUpTableG0[r][t]) / dtau;
		}

		void BuildHoppingMatrix()
		{
			for (uint_t i = 0; i < lattice.Sites(); ++i)
			{
				for (uint_t j = 0; j < lattice.Sites(); ++j)
				{
					if (lattice.IsNeighbor(i, j))
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
			lattice.Resize(l, rng);
			hoppingMatrix.resize(lattice.Sites(), lattice.Sites());
			lookUpTableG0.AllocateTable(lattice.MaxDistance() + 1, nTimeBins + 1);
			lookUpTableDtG0.AllocateTable(lattice.MaxDistance() + 1, nTimeBins);
			nhoodDist = std::min({uint_t(100), lattice.MaxDistance()});
		}
		
		void SetTemperature(value_t T)
		{
			beta = 1.0 / T;
			dtau = beta / static_cast<value_t>(nTimeBins);
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
		Geometry lattice;
		UpdateHandler_t updateHandler;
		value_t beta;
		value_t V;
		value_t t;
		value_t infinTau;
		value_t zeta2;
		value_t zeta4;
		uint_t nTimeBins;
		uint_t maxWorms = 4;
		LookUpTable<value_t, uint_t, 2> lookUpTableG0;
		LookUpTable<value_t, uint_t, 2> lookUpTableDtG0;
		value_t dtau;
		StateType state = StateType::Z;
		matrix_t hoppingMatrix;
		matrix_t hopDiag;
		matrix_t hopEV;
		matrix_t hopEVT;
		Eigen::SelfAdjointEigenSolver<matrix_t> evSolver;
		//Eigen::FullPivHouseholderQR<matrix_t> invSolver;
		Eigen::FullPivLU<matrix_t> invSolver;
		uint_t nhoodDist;
		std::list<std::string> updateList;
};