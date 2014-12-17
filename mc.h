#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <utility>
#include <cmath>
#include <numeric>
#include <string>
#include <chrono>
#include <cstdint>
#include <map>
#include "ConfigSpace.h"
#include "Random.h"
#include "HexagonalHoneycomb.h"
#include "RhombicHoneycomb.h"
//#define EIGEN_USE_MKL_ALL
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "measurements.h"
#include "random.h"
#include "parser.h"
#include "types.h"

using duration_t = std::chrono::seconds;

class make_string
{
public:
	template <typename T>
	make_string& operator<<(const T& val)
	{
		buffer << val;
		return *this;
	}
	operator std::string() const
	{
		return buffer.str();
	}
private:
	std::ostringstream buffer;
};

class mc
{
	public:
		using uint_t = std::uint_fast32_t;
		using int_t = std::int_fast32_t;
		using value_t = double;
		using matrix_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;
		using ConfigSpace_t = ConfigSpace<HexagonalHoneycomb<Random, uint_t, int_t>, Random, value_t, matrix_t>;

		mc(const std::string& dir);
		~mc();

		void random_write(odump& d);
		void seed_write(const std::string& fn);
		void random_read(idump& d);
		void init();
		void write(const std::string& dir);
		bool read(const std::string& dir);
		void write_output(const std::string& dir);
		bool is_thermalized();
		measurements measure;
		
	public:
		void do_update();
		void do_measurement();
		void BuildUpdateWeightMatrix();
		void PrintAcceptanceMatrix();
		void SelfBalance();
		void FinalizeSimulation();
		
	private:
		uint_t& GetWithDef(std::map<uint_t, uint_t>& map, uint_t key, uint_t defval)
		{
			auto it = map.find( key );
			if (it == map.end())
				map[key] = 0;
			return map[key];
		}
		
	private:
		Random rng;
		ConfigSpace_t configSpace;
		uint_t nThermalize;
		uint_t nMeasurements;
		uint_t nRebuild;
		uint_t nThermStep;
		uint_t nPrebins;
		int nUpdateType = 11;
		int nStateType = 3;
		matrix_t updateWeightMatrix = matrix_t(nUpdateType, nStateType);
		matrix_t acceptedUpdates = matrix_t(nUpdateType, nStateType);
		matrix_t proposedUpdates = matrix_t(nUpdateType, nStateType);
		matrix_t accRateUpdates = matrix_t(nUpdateType, nStateType);
		parser param;
		uint_t sweep = 0;
		uint_t rebuildCnt = 0;
		std::vector< value_t > corrVector;
		std::map<uint_t, uint_t> exporderHistZ;
		std::map<uint_t, uint_t> exporderHistW2;
		double* evalableParameters;
};
