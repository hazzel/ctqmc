#include "mc.h"
#include <sstream>
#include <fstream>
#include <string>
#include <limits>
#include <functional>
#include <omp.h>
#include <gperftools/profiler.h>
#include "glob.h"

std::vector<string> GlobFile(const string& pattern)
{
	glob_t glob_result;
	glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
	std::vector<string> files;
	for(unsigned int i=0;i<glob_result.gl_pathc;++i)
	{
		files.push_back(string(glob_result.gl_pathv[i]));
	}
	globfree(&glob_result);
	return files;
}

void M2Function(double& out, std::vector< std::valarray<double>* >& o, double* p)
{
	double z=(*o[0])[0];
	double w2=(*o[1])[0];
	double w4=(*o[2])[0];
	
	out = (w2 / z) / (p[0] * p[1] * p[3] * p[3]);
}

void M4Function(double& out, std::vector< std::valarray<double>* >& o, double* p)
{
	double z=(*o[0])[0];
	double w2=(*o[1])[0];
	double w4=(*o[2])[0];
	
	out = (w4 / z) / (p[0] * p[2] * p[3] * p[3] * p[3] * p[3]);
}

void BinderRatioFunction(double& out, std::vector< std::valarray<double>* >& o, double* p)
{
	double z=(*o[0])[0];
	double w2=(*o[1])[0];
	double w4=(*o[2])[0];
	
	out = (p[0] * p[1] * p[1] / p[2]) * (w4 * z / (w2 * w2));
}

void ChiFunction(double& out, std::vector< std::valarray<double>* >& o, double* p)
{
	double z=(*o[0])[0];
	double w2=(*o[1])[0];
	double w4=(*o[2])[0];
	
	out = (w2 / z) / (p[0] * p[1] * p[3] * p[3]);
}

void AvgExporderFunction(double& out, std::vector< std::valarray<double>* >& o)
{
	double k=(*o[0])[0];
	double z=(*o[1])[0];
	
	out = k / z;
}

void CorrFunction(std::valarray<double>& out, std::vector< std::valarray<double>* >& o, double* p)
{
	std::valarray<double>* corr = o[0];
	double z=(*o[1])[0];
	out.resize(corr->size());
	for (int i = 0; i < corr->size(); ++i)
		out[i] = (*corr)[i] / (p[0] * p[1] * z);
}

CLASSNAME::CLASSNAME(const std::string& dir)
	: rng(Random()), configSpace(rng)
{
	nUpdateType = 16;
	nStateType = 3;
	sweep = 0;
	rebuildCnt = 0;
	isInitialized = false;
	nOptimizationSteps = 0;
	annealing = false;
	startT = 5.0;
	updateWeightMatrix.resize(nUpdateType, nStateType);
	proposeProbabilityMatrix.resize(nUpdateType, nStateType);
	acceptedUpdates.resize(nUpdateType, nStateType);
	proposedUpdates.resize(nUpdateType, nStateType);

	std::cout.precision(15);
	std::cout << std::fixed;
	
	param.read_file(dir);
	value_t T = param.value_or_default<value_t>("T", 1.);
	L = param.value_or_default<uint_t>("L", 6);
	pt_spacing = param.value_or_default<uint_t>("GLOBAL_UPDATE_SPACING", 100);
	label = 0;
	
	#ifdef MCL_PT
		pt_var = param.return_vector<value_t>("@V");
		for (uint_t i=0; i < pt_var.size(); ++i)
		{
			measurements newmeasure;
			measure.push_back(newmeasure);
			exporderHistZ.push_back(std::map<uint_t, uint_t>());
			exporderHistW2.push_back(std::map<uint_t, uint_t>());
			exporderHistW4.push_back(std::map<uint_t, uint_t>());
		}
	#endif
	
	//configSpace.fileIO = param.value_or_default<uint_t>("FILEIO", 0);
	configSpace.fileIO = false;
	configSpace.nTimeBins = param.value_or_default<uint_t>("TIMEBINS", 50000);
	configSpace.t = param.value_or_default<value_t>("t0", 1.0);
	#ifndef MCL_PT
		configSpace.V = param.value_or_default<value_t>("V", 1.4);
	#endif
	finalT = T;
	//startT = T;
	if (annealing)
		configSpace.SetTemperature(startT);
	else
		configSpace.SetTemperature(T);
	std::string geometry = param.value_or_default<std::string>("GEOMETRY", "hex");
	path = dir.substr(0, dir.substr(0, dir.rfind('/')).rfind('/') + 1);
	
	if (geometry == "hex")
		std::cout << "error: need to reimplement hex" << std::endl;
	else if (geometry == "rhom")
	{
		std::string geo_file = path + "geometry/rhom-L" + ToString(L);
		configSpace.lattice = new Rhom_t(geo_file, true);
	}
	//Set up geometry
	configSpace.ResizeGeometry(L);
	configSpace.BuildHoppingMatrix();

	//Build G0 look up table
	//std::string g0_file = path + "g0lookup/" + geometry + "-B" + ToString(configSpace.nTimeBins / 1000) + "-L" + ToString(L) + "-T" + ToString(T);
	configSpace.BuildG0LookUpTable();
	configSpace.updateHandler.Init();
	
	int_t nhd = param.value_or_default<int_t>("NHOODDIST", 1);
	configSpace.nhoodDist = std::min(nhd, configSpace.lattice->MaxDistance());
	configSpace.zeta2 = param.value_or_default<value_t>("zeta2", 1.0);
	value_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
	//configSpace.zeta2 /= m / T;
	configSpace.zeta2 /= std::pow(L, 4.) / T;
	configSpace.zeta4 = param.value_or_default<value_t>("zeta4", 1.0);
	//configSpace.zeta4 /= m * m * m / T;
	configSpace.zeta4 /= std::pow(L, 8.)  / T;

	nOptimizationSteps = param.value_or_default<value_t>("ZETA_OPTIMIZATION", 0.0);
	nOptimizationTherm = param.value_or_default<value_t>("ZETA_THERM", 10000.0);

	nThermalize = param.value_or_default<uint_t>("THERMALIZATION", 10000);
	nMeasurements = param.value_or_default<uint_t>("SWEEPS", 10000);
	nPrebins = param.value_or_default<uint_t>("PREBINS", 100);
	nRebuild = param.value_or_default<uint_t>("REBUILD", 100000);
	nThermStep = param.value_or_default<uint_t>("THERM_STEP", 100);
	if (param.defined("SEED"))
		rng.NewRng(param.value_of<uint_t>("SEED"));
	
	evalableParameters = new double[4];
	evalableParameters[0] = configSpace.beta;
	evalableParameters[1] = configSpace.zeta2;
	evalableParameters[2] = configSpace.zeta4;
	evalableParameters[3] = configSpace.lattice->Sites();
	corrVector.resize(configSpace.lattice->MaxDistance() + 1, 0.0);
	
	bool do_therm = param.value_or_default<uint_t>("THERMALIZE", 0);
	//Read thermalized state if it exists
	therm_name_depr = "therm-L" + ToString(L) + "-V" + ToString(configSpace.V) + "-T" + ToString(T) + "-" + geometry;
	therm_name = therm_name_depr + "-N" + ToString(nThermalize);
	therm_path = path + "thermalization/";
	std::vector<std::string> files = GlobFile(therm_path + therm_name + "/*.txt");
	if (FileExists(therm_path + therm_name_depr) && (!DirExists(therm_path + therm_name)))
		files.push_back(therm_path + therm_name_depr);
	if (sweep == 0 && (!do_therm) && files.size() > 0 && nOptimizationSteps == 0)
	{
		read_state(files[rng() * files.size()]);
		sweep = nThermalize;
		std::cout << "Thermalization...Done" << std::endl;
	}
	
	BuildUpdateWeightMatrix();
	//ProfilerStart("gperf/mc.prof");
}

CLASSNAME::~CLASSNAME()
{
	delete[] evalableParameters;
	//ProfilerStop();
}

void CLASSNAME::random_write(odump& d)
{
	rng.RngHandle()->write(d);
}
void CLASSNAME::seed_write(const std::string& fn)
{
	std::ofstream s;
	s.open(fn.c_str());
	s << rng.Seed() << std::endl;
	s.close();
}
void CLASSNAME::random_read(idump& d)
{
	rng.NewRng();
	rng.RngHandle()->read(d);
}
void CLASSNAME::init()
{
	#ifdef MCL_PT
		for (uint_t i=0; i< pt_var.size(); ++i)
		{
			measure[i].add_observable("k", nPrebins);
			measure[i].add_observable("<w>", nPrebins);
			measure[i].add_observable("deltaZ", nPrebins);
			measure[i].add_observable("deltaW2", nPrebins);
			measure[i].add_observable("deltaW4", nPrebins);
			measure[i].add_observable("deltaChi", nPrebins);
			measure[i].add_observable("avgInvGError");
			measure[i].add_observable("condition");
			measure[i].add_observable("weight", nPrebins);
			measure[i].add_vectorobservable("Corr", configSpace.lattice->MaxDistance() + 1, nPrebins);
		}
	#else
		measure.add_observable("k", nPrebins);
		measure.add_observable("<w>", nPrebins);
		measure.add_observable("deltaZ", nPrebins);
		measure.add_observable("deltaW2", nPrebins);
		measure.add_observable("deltaW4", nPrebins);
		measure.add_observable("deltaChi", nPrebins);
		measure.add_observable("avgInvGError");
		measure.add_observable("condition");
		measure.add_observable("weight", nPrebins);
		measure.add_vectorobservable("Corr", configSpace.lattice->MaxDistance() + 1, nPrebins);
	#endif
}
void CLASSNAME::write(const std::string& dir)
{
	odump d(dir+"dump");
	random_write(d);
	d.write(sweep);
	//std::cout << "write sweep: " << sweep << std::endl;
	d.write(rebuildCnt);
	for (uint_t i = 0; i < nUpdateType; ++i)
	{
		for (uint_t j = 0; j < nStateType; ++j)
		{
			d.write(updateWeightMatrix(i, j));
			d.write(acceptedUpdates(i, j));
			d.write(proposedUpdates(i, j));
		}
	}
	configSpace.Serialize(d);
	d.close();
	seed_write(dir+"seed");
	std::ofstream f; f.open((dir + "bins").c_str());
	f << ( sweep > nThermalize ? sweep-nThermalize : 0 ) << std::endl;
	f.close();

	std::ofstream ostream;
	#ifdef MCL_PT
		ostream.open(dir+"exporderhist.para"+std::to_string(myrep+1)+".txt");
		for (uint_t i = 0; i < std::max(exporderHistZ[myrep].size(), exporderHistW2[myrep].size()); ++i)
			ostream << i << " " << GetWithDef(exporderHistZ[myrep], i, 0) << " " << GetWithDef(exporderHistW2[myrep], i, 0) << std::endl;
		ostream.close();
	#else
		std::string ofile(dir+"exporderhist.txt");
		ostream.open(ofile.c_str());
		for (uint_t i = 0; i < std::max(exporderHistZ.size(), exporderHistW2.size()); ++i)
			ostream << i << " " << GetWithDef(exporderHistZ, i, 0) << " " << GetWithDef(exporderHistW2, i, 0) << std::endl;
		ostream.close();
	#endif
	
	ostream.open((dir+"probabilities.txt").c_str());
	PrintAcceptanceMatrix(ostream);
	ostream.close();
}
void CLASSNAME::write_state(const std::string& dir)
{
	odump d(dir);
	configSpace.Serialize(d);
	d.close();
}
bool CLASSNAME::read(const std::string& dir)
{
	//std::cout << "read " << dir <<  std::endl;
	idump d(dir+"dump");
	if (!d)
	{
		std::cout << "read fail" << std::endl;
		return false;
	}
	else
	{
		random_read(d);
		d.read(sweep);
		d.read(rebuildCnt);
		for (uint_t i = 0; i < nUpdateType; ++i)
		{
			for (uint_t j = 0; j < nStateType; ++j)
			{
				d.read(updateWeightMatrix(i, j));
				d.read(acceptedUpdates(i, j));
				d.read(proposedUpdates(i, j));
			}
		}
		configSpace.Serialize(d);
		d.close();
		//std::cout << "read sweep: " << sweep << " , pertorder:" << configSpace.updateHandler.GetVertexHandler().Vertices() << std::endl;
		return true;
	}
}
bool CLASSNAME::read_state(const std::string& dir)
{
	std::ifstream is(dir);
	configSpace.SerializeTxt(is);
	is.close();
	return true;
}

#ifdef MCL_PT
void CLASSNAME::write_output(const std::string& dir, int para)
{
	measure[para].add_evalable("M2","deltaZ","deltaW2","deltaW4", M2Function, evalableParameters);
	measure[para].add_evalable("M4","deltaZ","deltaW2","deltaW4", M4Function, evalableParameters);
	measure[para].add_evalable("BinderRatio","deltaZ","deltaW2","deltaW4", BinderRatioFunction, evalableParameters);
	measure[para].add_evalable("Chi","deltaZ","deltaChi","deltaW4", ChiFunction, evalableParameters);
	//measure[para].add_evalable("AvgExpOrder","k","deltaZ", AvgExporderFunction);
	measure[para].add_vectorevalable("Correlations","Corr","deltaZ", CorrFunction, evalableParameters);
	std::ofstream f;
	//f.precision(8);
	//f << std::fixed;
	f.open(dir.c_str());
	f << "PARAMETERS" << endl;
	param.get_all_with_one_from_specified_array("@V", para, f);
	measure[para].get_statistics(f);
}
#else
void CLASSNAME::write_output(const std::string& dir)
{
	measure.add_evalable("M2","deltaZ","deltaW2","deltaW4", M2Function, evalableParameters);
	measure.add_evalable("M4","deltaZ","deltaW2","deltaW4", M4Function, evalableParameters);
	measure.add_evalable("BinderRatio","deltaZ","deltaW2","deltaW4", BinderRatioFunction, evalableParameters);
	measure.add_evalable("Chi","deltaZ","deltaChi","deltaW4", ChiFunction, evalableParameters);
	//measure.add_evalable("AvgExpOrder","k","deltaZ", AvgExporderFunction);
	measure.add_vectorevalable("Correlations","Corr","deltaZ", CorrFunction, evalableParameters);
	std::ofstream f;
	//f.precision(8);
	//f << std::fixed;
	f.open(dir.c_str());
	f << "PARAMETERS" << endl;
	param.get_all(f);
	measure.get_statistics(f);
}
#endif

bool CLASSNAME::is_thermalized()
{
	return (nZetaOptimization >= nOptimizationSteps) && (sweep >= nThermalize);
}

void CLASSNAME::BuildUpdateWeightMatrix()
{
	
	/*
	//Z<->W2 LEI
	proposeProbabilityMatrix <<	2.05 / 10.0	,	1.30 / 10.0	,	1.65 / 10.0,
															2.05 / 10.0	,	1.30 / 10.0	,	1.65 / 10.0,
															1.5 / 10.0	,	1.0 / 10.0	,	1.0 / 10.0,
															1.5 / 10.0	,	1.0 / 10.0	,	1.0 / 10.0,
															0.4 / 10.0	,	0.25 / 10.0	,	0.5 / 10.0,
															0.4 / 10.0	,	0.25 / 10.0	,	0.5 / 10.0,
															0.4 / 10.0	,	0.25 / 10.0	,	0.25 / 10.0,
															0.4 / 10.0	,	0.25 / 10.0	,	0.25 / 10.0,
															1.3 / 10.0	,	0.0					,	0.0,
															0.0					,	1.8 / 10.0	,	0.0, 
															0.0 / 10.0	,	0.0					,	0.0,
															0.0					,	0.0					,	0.0 / 10.0,
															0.0					,	0.6 / 10.0	,	0.0,
															0.0					,	0.0					,	1.2 / 10.0,
															0.0					,	0.0 / 10.0	,	0.0 / 10.0,
															0.0					,	2.0 / 10.0	,	2.0 / 10.0;
			*/


	//ALL TRANSITIONS LEI
	proposeProbabilityMatrix <<	2.0 / 10.0	,	1.30 / 10.0	,	1.25 / 10.0,
															2.0 / 10.0	,	1.30 / 10.0	,	1.25 / 10.0,
															1.5 / 10.0	,	1.0 / 10.0	,	1.0 / 10.0,
															1.5 / 10.0	,	1.0 / 10.0	,	1.0 / 10.0,
															0.25 / 10.0	,	0.25 / 10.0	,	0.5 / 10.0,
															0.25 / 10.0	,	0.25 / 10.0	,	0.5 / 10.0,
															0.25 / 10.0	,	0.25 / 10.0	,	0.25 / 10.0,
															0.25 / 10.0	,	0.25 / 10.0	,	0.25 / 10.0,
															1.2 / 10.0	,	0.0					,	0.0,
															0.0					,	1.8 / 10.0	,	0.0, 
															0.8 / 10.0	,	0.0					,	0.0,
															0.0					,	0.0					,	1.2 / 10.0,
															0.0					,	0.6 / 10.0	,	0.0,
															0.0					,	0.0					,	1.2 / 10.0,
															0.0					,	0.0 / 10.0	,	0.0 / 10.0,
															0.0					,	2.0 / 10.0	,	1.6 / 10.0;


	
	/*
	//ALL TRANSITIONS LEI
	proposeProbabilityMatrix <<	1.0 / 10.0	,	1.0 / 10.0	,	1.0 / 10.0,
															1.0 / 10.0	,	1.0 / 10.0	,	1.0 / 10.0,
															0.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
															0.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
															0.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
															0.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
															0.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
															0.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
															4.0 / 10.0	,	0.0					,	0.0,
															0.0					,	3.0 / 10.0	,	0.0, 
															4.0 / 10.0	,	0.0					,	0.0,
															0.0					,	0.0					,	3.0 / 10.0,
															0.0					,	3.0 / 10.0	,	0.0,
															0.0					,	0.0					,	3.0 / 10.0,
															0.0					,	0.0 / 10.0	,	0.0 / 10.0,
															0.0					,	2.0 / 10.0	,	2.0 / 10.0;
	*/
	/*
		//ALL TRANSITIONS
	proposeProbabilityMatrix <<	1.05 / 10.0	,	0.8 / 10.0	,	0.7 / 10.0,
											3.05 / 10.0	,	2.4 / 10.0	,	2.3 / 10.0,
											0.5 / 10.0	,	0.5 / 10.0	,	0.5 / 10.0,
											1.5 / 10.0	,	1.5 / 10.0	,	1.5 / 10.0,
											0.5 / 10.0	,	0.25 / 10.0	,	0.25 / 10.0,
											1.5 / 10.0	,	0.25 / 10.0	,	0.25 / 10.0,
											0.25 / 10.0	,	0.25 / 10.0	,	0.25 / 10.0,
											0.75 / 10.0	,	0.25 / 10.0	,	0.25 / 10.0,
											0.8 / 10.0	,	0.0			,	0.0,
											0.0			,	1.2 / 10.0	,	0.0, 
											0.1 / 10.0	,	0.0			,	0.0,
											0.0			,	0.0			,	1.2 / 10.0,
											0.0			,	0.6 / 10.0	,	0.0,
											0.0			,	0.0			,	1.2 / 10.0,
											0.0			,	0.0 / 10.0	,	0.0 / 10.0,
											0.0			,	2.0 / 10.0	,	1.6 / 10.0;
	*/
	
	/*
	//ALL TRANSITIONS
	proposeProbabilityMatrix <<	1.0 / 10.0,	1. / 10.0,	1. / 10.0,
															1.0 / 10.0,	1. / 10.0,	1. / 10.0,
															0.25 / 10.0,	0.25 / 10.0,	0.25 / 10.0,
															0.25 / 10.0,	0.25 / 10.0,	0.25 / 10.0,
															0.25 / 10.0,	0.25 / 10.0,	0.25 / 10.0,
															0.25 / 10.0,	0.25 / 10.0,	0.25 / 10.0,
															0.25 / 10.0,	0.25 / 10.0,	0.25 / 10.0,
															0.25 / 10.0,	0.25 / 10.0,	0.25 / 10.0,
															3.25 / 10.0,	0.0				,	0.0,
															0.0				,	2.75 / 10.0,	0.0, 
															3.25 / 10.0,	0.0				,	0.0,
															0.0				,	0.0				,	2.75 / 10.0,
															0.0				,	2.75 / 10.0,	0.0,
															0.0				,	0.0				,	2.75 / 10.0,
															0.0				,	0.0 / 10.0,	0.0 / 10.0,
															0.0				,	1.0 / 10.0,	1.0 / 10.0;
*/
	
	/*
		//TEST
	proposeProbabilityMatrix <<	5.0 / 10.0	,	5.0 / 10.0	,	5.0 / 10.0,
											10.0 / 10.0	,	10.0 / 10.0	,	10.0 / 10.0,
											1.0 / 10.0	,	0.5 / 10.0	,	0.0 / 10.0,
											1.0 / 10.0	,	0.5 / 10.0	,	0.0 / 10.0,
											0.5 / 10.0	,	0.5 / 10.0	,	0.0 / 10.0,
											0.5 / 10.0	,	0.5 / 10.0	,	0.0 / 10.0,
											0.5 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
											0.5 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
											2.0 / 10.0	,	0.0			,	0.0,
											0.0			,	2.0 / 10.0	,	0.0, 
											2.0 / 10.0	,	0.0			,	0.0,
											0.0			,	0.0			,	2.0 / 10.0,
											0.0			,	2.0 / 10.0	,	0.0,
											0.0			,	0.0			,	2.0 / 10.0,
											0.0			,	0.0 / 10.0	,	0.0 / 10.0,
											0.0			,	2.0 / 10.0	,	0.0 / 10.0;
	*/
	
	updateWeightMatrix.setZero();
	for (uint_t i = 0; i < nUpdateType; ++i)
	{
		for (uint_t j = 0; j < nStateType; ++j)
		{
			if (i == 0)
				updateWeightMatrix(i, j) = proposeProbabilityMatrix(i, j);
			else
				updateWeightMatrix(i, j) = updateWeightMatrix(i - 1, j) + proposeProbabilityMatrix(i, j);
		}
	}
	acceptedUpdates.setZero();
	proposedUpdates.setZero();
}

void CLASSNAME::PrintAcceptanceMatrix(std::ostream& out)
{
	out << "Acceptance of updates:" << std::endl;
	for (uint_t i = 0; i < nUpdateType; ++i)
	{
		for (uint_t j = 0; j < nStateType; ++j)
		{
			if (proposedUpdates(i, j) == 0)
			{
				out << 0 << " ";
			}
			else
			{
				out << acceptedUpdates(i, j) / proposedUpdates(i, j) << " ";
			}
		}
		out << std::endl;
	}
	out << "Total updates accepted: " << acceptedUpdates.sum() << " (" << acceptedUpdates.sum() / proposedUpdates.sum() * 100.0 << " %)" << std::endl;
	out << "Total updates proposed: " << proposedUpdates.sum() << std::endl;
}

void CLASSNAME::do_update()
{
	if (sweep == 0 && nZetaOptimization == 0)
	{
		std::cout << std::endl;
		std::cout << "Thermalization" << std::endl;
	}

	for (uint_t i = 0; i < nThermStep; ++i)
	{
		value_t r = rng();
		StateType state = configSpace.State();
		
		if (r < updateWeightMatrix(AddVertex, state))
		{
			const int_t N = 1;
			value_t proposeRatio = proposeProbabilityMatrix(RemoveVertex, state) / proposeProbabilityMatrix(AddVertex, state);
			if (configSpace.AddRandomVertices<N>(proposeRatio, false))
			{
				acceptedUpdates(AddVertex, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(AddVertex, state) += 1.0;
		}
		else if (r < updateWeightMatrix(RemoveVertex, state))
		{
			const int_t N = 1;
			value_t proposeRatio = proposeProbabilityMatrix(AddVertex, state) / proposeProbabilityMatrix(RemoveVertex, state);
			if (configSpace.RemoveRandomVertices<N>(proposeRatio, false))
			{
				acceptedUpdates(RemoveVertex, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(RemoveVertex, state) += 1.0;
		}
		else if (r < updateWeightMatrix(Add2Vertices, state))
		{
			const int_t N = 2;
			value_t proposeRatio = proposeProbabilityMatrix(Remove2Vertices, state) / proposeProbabilityMatrix(Add2Vertices, state);
			if (configSpace.AddRandomVertices<N>(proposeRatio, false))
			{
				acceptedUpdates(Add2Vertices, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(Add2Vertices, state) += 1.0;
		}
		else if (r < updateWeightMatrix(Remove2Vertices, state))
		{
			const int_t N = 2;
			value_t proposeRatio = proposeProbabilityMatrix(Add2Vertices, state) / proposeProbabilityMatrix(Remove2Vertices, state);
			if (configSpace.RemoveRandomVertices<N>(proposeRatio, false))
			{
				acceptedUpdates(Remove2Vertices, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(Remove2Vertices, state) += 1.0;
		}
		else if (r < updateWeightMatrix(Add5Vertices, state))
		{
			const int_t N = 3;
			value_t proposeRatio = proposeProbabilityMatrix(Remove5Vertices, state) / proposeProbabilityMatrix(Add5Vertices, state);
			if (configSpace.AddRandomVertices<N>(proposeRatio, false))
			{
				acceptedUpdates(Add5Vertices, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(Add5Vertices, state) += 1.0;
		}
		else if (r < updateWeightMatrix(Remove5Vertices, state))
		{
			const int_t N = 3;
			value_t proposeRatio = proposeProbabilityMatrix(Add5Vertices, state) / proposeProbabilityMatrix(Remove5Vertices, state);
			if (configSpace.RemoveRandomVertices<N>(proposeRatio, false))
			{
				acceptedUpdates(Remove5Vertices, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(Remove5Vertices, state) += 1.0;
		}
		else if (r < updateWeightMatrix(Add8Vertices, state))
		{
			const int_t N = 4;
			value_t proposeRatio = proposeProbabilityMatrix(Remove8Vertices, state) / proposeProbabilityMatrix(Add8Vertices, state);
			if (configSpace.AddRandomVertices<N>(proposeRatio, false))
			{
				acceptedUpdates(Add8Vertices, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(Add8Vertices, state) += 1.0;
		}
		else if (r < updateWeightMatrix(Remove8Vertices, state))
		{
			const int_t N = 4;
			value_t proposeRatio = proposeProbabilityMatrix(Add8Vertices, state) / proposeProbabilityMatrix(Remove8Vertices, state);
			if (configSpace.RemoveRandomVertices<N>(proposeRatio, false))
			{
				acceptedUpdates(Remove8Vertices, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(Remove8Vertices, state) += 1.0;
		}
		else if (r < updateWeightMatrix(ZtoW2, state) && state == Z)
		{
			value_t proposeRatio = proposeProbabilityMatrix(W2toZ, W2) / proposeProbabilityMatrix(ZtoW2, Z);
			bool result;
			if (configSpace.rng() < 0.5)
				result = configSpace.AddRandomVertices<1>(proposeRatio, true);
			else
				result = configSpace.OpenUpdate<1>(proposeRatio);
			if (result)
			{
				acceptedUpdates(ZtoW2, state) += 1.0;
				configSpace.state = W2;
				++rebuildCnt;
			}
			proposedUpdates(ZtoW2, state) += 1.0;
		}
		else if (r < updateWeightMatrix(W2toZ, state) && state == W2)
		{
			value_t proposeRatio = proposeProbabilityMatrix(ZtoW2, Z) / proposeProbabilityMatrix(W2toZ, W2);
			bool result;
			if (configSpace.rng() < 0.5)
				result = configSpace.RemoveRandomVertices<1>(proposeRatio, true);
			else
				result = configSpace.CloseUpdate<1>(proposeRatio);
			if (result)
			{
				acceptedUpdates(W2toZ, state) += 1.0;
				configSpace.state = Z;
				++rebuildCnt;
			}
			proposedUpdates(W2toZ, state) += 1.0;
		}
		else if (r < updateWeightMatrix(ZtoW4, state) && state == Z)
		{
			value_t proposeRatio = proposeProbabilityMatrix(W4toZ, W4) / proposeProbabilityMatrix(ZtoW4, Z);
			if (configSpace.AddRandomVertices<2>(proposeRatio, true))
			{
				acceptedUpdates(ZtoW4, state) += 1.0;
				configSpace.state = W4;
				++rebuildCnt;
			}
			proposedUpdates(ZtoW4, state) += 1.0;
		}
		else if (r < updateWeightMatrix(W4toZ, state) && state == W4)
		{
			value_t proposeRatio = proposeProbabilityMatrix(ZtoW4, Z) / proposeProbabilityMatrix(W4toZ, W4);
			if (configSpace.RemoveRandomVertices<2>(proposeRatio, true))
			{
				acceptedUpdates(W4toZ, state) += 1.0;
				configSpace.state = Z;
				++rebuildCnt;
			}
			proposedUpdates(W4toZ, state) += 1.0;
		}
		else if (r < updateWeightMatrix(W2toW4, state) && state == W2)
		{
			value_t proposeRatio = proposeProbabilityMatrix(W4toW2, W4) / proposeProbabilityMatrix(W2toW4, W2);
			if (configSpace.AddRandomVertices<1>(proposeRatio, true))
			{
				acceptedUpdates(W2toW4, state) += 1.0;
				configSpace.state = W4;
				++rebuildCnt;
			}
			proposedUpdates(W2toW4, state) += 1.0;
		}
		else if (r < updateWeightMatrix(W4toW2, state) && state == W4)
		{
			value_t proposeRatio = proposeProbabilityMatrix(W2toW4, W2) / proposeProbabilityMatrix(W4toW2, W4);
			if (configSpace.RemoveRandomVertices<1>(proposeRatio, true))
			{
				acceptedUpdates(W4toW2, state) += 1.0;
				configSpace.state = W2;
				++rebuildCnt;
			}
			proposedUpdates(W4toW2, state) += 1.0;
		}
		else if (r < updateWeightMatrix(replaceWorm, state) && state != Z)
		{
			bool result;
			if (state == W2)
			{
				result = configSpace.ReplaceWorm<1>();
			}
			else if (state == W4)
			{
				result = false;
			}
			if (result)
			{
				acceptedUpdates(replaceWorm, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(replaceWorm, state) += 1.0;
		}
		else if (r < updateWeightMatrix(shiftWorm, state) && state != Z)
		{
			bool result;
			if (state == W2)
			{
				result = configSpace.ShiftWorm<1>();
				//if (!result)
				//	result = configSpace.ReplaceWorm<1>();
			}
			else if (state == W4)
				result = configSpace.ShiftWorm<2>();
			if (result)
			{
				acceptedUpdates(shiftWorm, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(shiftWorm, state) += 1.0;
		}

		if (rebuildCnt == nRebuild)
		{
			double err = configSpace.updateHandler.IsStableInverse();
			#ifdef MCL_PT
				measure[myrep].add("avgInvGError", err);
			#else
				measure.add("avgInvGError", err);
				measure.add("condition", 0.0);
			#endif
			rebuildCnt = 0;
		}
	}
	++sweep;
	if (sweep % 50000 == 0)
		std::cout << "sweep : " << sweep << " , ( " << static_cast<double>(sweep - nThermalize) * 100. / static_cast<double>(nMeasurements) << " % )" << std::endl;

	if (nZetaOptimization < nOptimizationSteps)
	{
		OptimizeZeta();
	}
	else
	{
		if (sweep == nThermalize)
		{
			std::cout << "Done" << std::endl;
			ClearExpOrderHist();
			if (!DirExists(therm_path + therm_name))
				if (mkdir((therm_path + therm_name).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
					std::cout << "Could not create directory." << std::endl;
			std::vector<std::string> files = GlobFile(therm_path + therm_name + "/*");
			write_state(therm_path + therm_name + "/" + therm_name + "-" + ToString(files.size()));
		}
	}
	if (annealing && !is_thermalized() && (sweep % (nThermalize / 25)) == 0)
	{
		ThermalizationTemp();
	}
}

void CLASSNAME::MeasureExpOrder()
{
	switch (configSpace.State())
	{
		case Z:
			#ifdef MCL_PT
				GetWithDef(exporderHistZ[myrep], configSpace.updateHandler.GetVertexHandler().Vertices(), 0) += 1;
			#else
				GetWithDef(exporderHistZ, configSpace.updateHandler.GetVertexHandler().Vertices(), 0) += 1;
			#endif
			break;
		case W2:
			#ifdef MCL_PT
				GetWithDef(exporderHistW2[myrep], configSpace.updateHandler.GetVertexHandler().Vertices(), 0) += 1;
			#else
				GetWithDef(exporderHistW2, configSpace.updateHandler.GetVertexHandler().Vertices(), 0) += 1;
			#endif
			break;
		case W4:
			#ifdef MCL_PT
				GetWithDef(exporderHistW4[myrep], configSpace.updateHandler.GetVertexHandler().Vertices(), 0) += 1;
			#else
				GetWithDef(exporderHistW4, configSpace.updateHandler.GetVertexHandler().Vertices(), 0) += 1;
			#endif
			break;
	}
}

void CLASSNAME::ClearExpOrderHist()
{
	#ifdef MCL_PT
		for (uint_t i = 0; i < pt_var.size(); ++i)
		{
			exporderHistZ[i].clear();
			exporderHistW2[i].clear();
			exporderHistW4[i].clear();
		}
	#else
		exporderHistZ.clear();
		exporderHistW2.clear();
		exporderHistW4.clear();
	#endif
}

void CLASSNAME::OptimizeZeta()
{
	if (sweep < nOptimizationTherm)
	{
		switch (configSpace.State())
		{
			case Z:
				therm.State[Z] = therm.State[Z] * therm.N / (therm.N + 1.0) + 1.0 / (therm.N + 1.0);
				therm.State[W2] = therm.State[W2] * therm.N / (therm.N + 1.0);
				therm.State[W4] = therm.State[W4] * therm.N / (therm.N + 1.0);
				break;
			case W2:
				therm.State[Z] = therm.State[Z] * therm.N / (therm.N + 1.0);
				therm.State[W2] = therm.State[W2] * therm.N / (therm.N + 1.0) + 1.0 / (therm.N + 1.0);
				therm.State[W4] = therm.State[W4] * therm.N / (therm.N + 1.0);
				break;
			case W4:
				therm.State[Z] = therm.State[Z] * therm.N / (therm.N + 1.0);
				therm.State[W2] = therm.State[W2] * therm.N / (therm.N + 1.0);
				therm.State[W4] = therm.State[W4] * therm.N / (therm.N + 1.0) + 1.0 / (therm.N + 1.0);
				break;
		}
		therm.N += 1.0;
	}
	else if (sweep == nOptimizationTherm)
	{
		configSpace.zeta2 += (1./3. - therm.State[W2]) * configSpace.zeta2;
		configSpace.zeta4 += (1./3. - therm.State[W4]) * configSpace.zeta4;
		
		configSpace.zeta2 *= (2./3. + therm.State[Z]);
		configSpace.zeta2 *= (2./3. + therm.State[W4]) / (2./3. + therm.State[W2]);
		
		configSpace.zeta4 *= (2./3. + therm.State[Z]);
		configSpace.zeta4 *= (2./3. + therm.State[W2]) / (2./3. + therm.State[W4]);
		
		evalableParameters[1] = configSpace.zeta2;
		evalableParameters[2] = configSpace.zeta4;
		prevZeta.push_back(std::make_pair(configSpace.zeta2, configSpace.zeta4));
		value_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
		therm.Reset();
		++nZetaOptimization;
		sweep = 0;
		std::cout << nZetaOptimization << std::endl;

		if (nZetaOptimization > 0)
		{
			value_t zeta2mean = 0.0, zeta4mean = 0.0;
			for (uint_t i = prevZeta.size() * 3 / 4; i < prevZeta.size(); ++i)
			{
				zeta2mean += prevZeta[i].first / static_cast<value_t>(prevZeta.size() / 4);
				zeta4mean += prevZeta[i].second / static_cast<value_t>(prevZeta.size() / 4);
			}
			if (nZetaOptimization == nOptimizationSteps)
			{
				sweep = nOptimizationSteps * nOptimizationTherm;
				configSpace.zeta2 = zeta2mean;
				configSpace.zeta4 = zeta4mean;
				evalableParameters[1] = configSpace.zeta2;
				evalableParameters[2] = configSpace.zeta4;
			}
			std::cout << "Zeta2(T=" << 1./configSpace.beta << ",V=" << configSpace.V << ") = " << configSpace.zeta2 * std::pow(configSpace.L, 4.) * configSpace.beta << std::endl;
			std::cout << "Zeta4(T=" << 1./configSpace.beta << ",V=" << configSpace.V << ") = " << configSpace.zeta4 * std::pow(configSpace.L, 8.) * configSpace.beta << std::endl;
		}
	}
}

void CLASSNAME::do_measurement()
{
	if ((sweep - nThermalize + 1) % (nMeasurements / 3) == 0)
	{
		std::cout << "Vertices: " << configSpace.updateHandler.GetVertexHandler().Vertices() << std::endl;
		std::cout << "Worms: " << configSpace.updateHandler.GetVertexHandler().Worms() << std::endl;
	}
	
	#ifdef MCL_PT
		measure[myrep].add("<w>", configSpace.updateHandler.GetVertexHandler().Worms());
		measure[myrep].add("weight", configSpace.updateHandler.GetWeight());
		measure[myrep].add("k", configSpace.updateHandler.GetVertexHandler().Vertices());
	#else
		measure.add("<w>", configSpace.updateHandler.GetVertexHandler().Worms());
		measure.add("weight", configSpace.updateHandler.GetWeight());
		measure.add("k", configSpace.updateHandler.GetVertexHandler().Vertices());
	#endif
	uint_t R = 0;
	value_t sign, c;
	std::fill(corrVector.begin(), corrVector.end(), 0.0);
	switch (configSpace.State())
	{
		case Z:
			#ifdef MCL_PT
				measure[myrep].add("deltaZ", 1.0);
				measure[myrep].add("deltaW2", 0.0);
				measure[myrep].add("deltaW4", 0.0);
				measure[myrep].add("deltaChi", 0.0);
			#else
				measure.add("deltaZ", 1.0);
				measure.add("deltaW2", 0.0);
				measure.add("deltaW4", 0.0);
				measure.add("deltaChi", 0.0);
			#endif
			break;

		case W2:
			sign = configSpace.updateHandler.GetVertexHandler().WormParity();
			#ifdef MCL_PT
				measure[myrep].add("deltaZ", 0.0);
				measure[myrep].add("deltaW2", 1.0);
				measure[myrep].add("deltaW4", 0.0);
				measure[myrep].add("deltaChi", sign);
			#else
				measure.add("deltaZ", 0.0);
				measure.add("deltaW2", 1.0);
				measure.add("deltaW4", 0.0);
				measure.add("deltaChi", sign);
			#endif

			R = configSpace.updateHandler.GetVertexHandler().WormDistance();
			corrVector[R] = sign / configSpace.lattice->DistanceHistogram(R);
			break;
		case W4:
			#ifdef MCL_PT
				measure[myrep].add("deltaZ", 0.0);
				measure[myrep].add("deltaW2", 0.0);
				measure[myrep].add("deltaW4", 1.0);
				measure[myrep].add("deltaChi", 0.0);
			#else
				measure.add("deltaZ", 0.0);
				measure.add("deltaW2", 0.0);
				measure.add("deltaW4", 1.0);
				measure.add("deltaChi", 0.0);
			#endif
			break;
	}
	#ifdef MCL_PT
		measure[myrep].add("Corr", corrVector);
	#else
		measure.add("Corr", corrVector);
	#endif
}

#ifdef MCL_PT
void CLASSNAME::change_to (int i)
{
	change_parameter(i);
}

void CLASSNAME::change_parameter(int i)
{
	myrep=i;
	configSpace.V = pt_var[myrep];
	if (myrep == 0)
		label = 1;
	else if (myrep == (int) (pt_var.size() - 1))
		label = -1;
}

bool CLASSNAME::request_global_update()
{
	bool result = (sweep && (sweep%pt_spacing==0)) ? true : false;
	return result;
}

double CLASSNAME::get_weight(int f)
{
	uint_t k = configSpace.updateHandler.GetVertexHandler().Vertices();
	return k * std::log(pt_var[f] / configSpace.V);
}


int CLASSNAME::get_label()
{
	return label;
}
#endif
