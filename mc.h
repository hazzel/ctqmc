#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <cmath>
#include <numeric>
#include <string>
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

#ifdef MCL_PT
	#define CLASSNAME mc_pt
#else
	#define CLASSNAME mc
#endif

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

struct Measure
{
	double N;
	double State[3];

	void Reset()
	{
		N = 0.0;
		State[0] = 0.0;
		State[1] = 0.0;
		State[2] = 0.0;
	}

	double Sigma()
	{
		return std::sqrt(std::pow((State[0] - 1./3.), 2.0) + std::pow((State[1] - 1./3.), 2.0) + std::pow((State[2] - 1./3.), 2.0));
	}
};

class CLASSNAME
{
	public:
		using uint_t = std::int_fast32_t;
		using int_t = std::int_fast32_t;
		using value_t = double;
		using matrix_t = Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>;
		using Hex_t = HexagonalHoneycomb<Random, int_t>;
		using Rhom_t = RhombicHoneycomb<Random, int_t>;
		using ConfigSpace_t = ConfigSpace<GeometryBase<Random, int_t>, Random, value_t, matrix_t>;

		CLASSNAME(const std::string& dir);
		~CLASSNAME();

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
		void PrintAcceptanceMatrix(std::ostream& out);
		void SelfBalance();
		void FinalizeSimulation();
		void OptimizeZeta();
		
	private:
		uint_t& GetWithDef(std::map<uint_t, uint_t>& map, uint_t key, uint_t defval)
		{
			auto it = map.find( key );
			if (it == map.end())
				map[key] = defval;
			return map[key];
		}
		
		value_t ThermalizationTemp()
		{
			value_t T = finalT + (nThermalize - sweep) / nThermalize * (startT - finalT);
			configSpace.SetTemperature(T);
			evalableParameters[0] = configSpace.beta;
			configSpace.BuildG0LookUpTable("");
		}
		
	private:
		Random rng;
		ConfigSpace_t configSpace;
		uint_t nThermalize;
		uint_t nMeasurements;
		uint_t nRebuild;
		uint_t nThermStep;
		uint_t nPrebins;
		int nUpdateType = 15;
		int nStateType = 3;
		matrix_t updateWeightMatrix = matrix_t(nUpdateType, nStateType);
		matrix_t acceptedUpdates = matrix_t(nUpdateType, nStateType);
		matrix_t proposedUpdates = matrix_t(nUpdateType, nStateType);
		parser param;
		uint_t sweep = 0;
		uint_t rebuildCnt = 0;
		std::vector< value_t > corrVector;
		std::map<uint_t, uint_t> exporderHistZ;
		std::map<uint_t, uint_t> exporderHistW2;
		std::map<uint_t, uint_t> exporderHistW4;
		double* evalableParameters;
		uint_t L;
		bool isInitialized = false;
		std::string path;
		Measure therm;
		std::map< value_t, std::pair<value_t, value_t> > zetaOptimization;
		uint_t nZetaOptimization = 0;
		uint_t nOptimizationSteps;
		uint_t nOptimizationTherm;
		bool annealing = true;
		value_t startT = 2.5;
		value_t finalT;
};
