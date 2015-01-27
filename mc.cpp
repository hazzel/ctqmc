#include "mc.h"
#include <sstream>
#include <fstream>
#include <string>
#include <limits>
#include <functional>

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

mc::mc(const std::string& dir)
	: rng(Random()), configSpace(rng)
{
	//fpu_fix_start(&old_cw);
	
	std::cout.precision(15);
	std::cout << std::fixed;
	
	param.read_file(dir);
	value_t T = param.value_or_default<value_t>("T", 1.);
	L = param.value_or_default<uint_t>("L", 4);
	
	configSpace.nTimeBins = param.value_or_default<uint_t>("TIMEBINS", 50000);
	configSpace.t = param.value_or_default<value_t>("t0", 1.0);
	configSpace.V = param.value_or_default<value_t>("V", 1.4);
	configSpace.SetTemperature(T);
	std::string geometry = param.value_or_default<std::string>("GEOMETRY", "hex");
	path = dir.substr(0, dir.rfind('/')+1);
	if (geometry == "hex")
		configSpace.lattice = new Hex_t();
	else if (geometry == "rhom")
	{
		std::string filename = path + "geometry/rhom-L" + std::to_string(L);
		configSpace.lattice = new Rhom_t(filename);
	}
	std::cout << "Set up geometry...";
	std::cout.flush();
	configSpace.ResizeGeometry(L);
	configSpace.BuildHoppingMatrix();

	std::cout << "Done." << std::endl;
	
	configSpace.zeta2 = param.value_or_default<value_t>("zeta2", 1.0);
	value_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
	configSpace.zeta2 /= m / T;
	configSpace.zeta4 = param.value_or_default<value_t>("zeta4", 1.0);
	configSpace.zeta4 /= m * m * m / T;
	
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
	
	BuildUpdateWeightMatrix();
}

mc::~mc()
{
	delete[] evalableParameters;
	//fpu_fix_end(&old_cw);
}

void mc::random_write(odump& d)
{
	rng.RngHandle()->write(d);
}
void mc::seed_write(const std::string& fn)
{
	std::ofstream s;
	s.open(fn.c_str());
	s << rng.Seed() << std::endl;
	s.close();
}
void mc::random_read(idump& d)
{
	rng.NewRng();
	rng.RngHandle()->read(d);
}
void mc::init()
{
	measure.add_observable("k", nPrebins);
	measure.add_observable("<w>", nPrebins);
	measure.add_observable("deltaZ", nPrebins);
	measure.add_observable("deltaW2", nPrebins);
	measure.add_observable("deltaW4", nPrebins);
	measure.add_observable("avgInvGError", nPrebins);
	measure.add_observable("relInvGError", nPrebins);
	measure.add_observable("condition", nPrebins);
	measure.add_vectorobservable("Corr", configSpace.lattice->MaxDistance() + 1, nPrebins);
}
void mc::write(const std::string& dir)
{
	odump d(dir+"dump");
	random_write(d);
	d.write(sweep);
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
	std::string ofile(dir+"exporderhist.txt");
	std::ofstream ostream;
	ostream.open(ofile.c_str());
	for (uint_t i = 0; i < std::max(exporderHistZ.size(), exporderHistW2.size()); ++i)
		ostream << i << " " << GetWithDef(exporderHistZ, i, 0) << " " << GetWithDef(exporderHistW2, i, 0) << std::endl;
	ostream.close();
	ostream.open(dir+"probabilities.txt");
	PrintAcceptanceMatrix(ostream);
	ostream.close();
}
bool mc::read(const std::string& dir)
{
	idump d(dir+"dump");
	if (!d) 
		return false;
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
		return true;
	}
}

void mc::write_output(const std::string& dir)
{
	measure.add_evalable("M2","deltaZ","deltaW2","deltaW4", M2Function, evalableParameters);
	measure.add_evalable("M4","deltaZ","deltaW2","deltaW4", M4Function, evalableParameters);
	measure.add_evalable("BinderRatio","deltaZ","deltaW2","deltaW4", BinderRatioFunction, evalableParameters);
	measure.add_evalable("AvgExpOrder","k","deltaZ", AvgExporderFunction);
	measure.add_vectorevalable("Correlations","Corr","deltaZ", CorrFunction, evalableParameters);
	std::ofstream f;
	//f.precision(8);
	//f << std::fixed;
	f.open(dir.c_str());
	f << "PARAMETERS" << endl;
	param.get_all(f);
	measure.get_statistics(f);
}

bool mc::is_thermalized()
{
	return sweep >= nThermalize;
}

void mc::BuildUpdateWeightMatrix()
{
/*
	//ALL TRANSITIONS
	updateWeightMatrix <<	1.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
									2.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
									3.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
									4.0 / 10.0	,	0.0 / 10.0	,	0.0 / 10.0,
									7.0 / 10.0	,	0.0			,	0.0,
									0.0			,	3.0 / 10.0	,	0.0, 
									10.0 / 10.0	,	0.0			,	0.0,
									0.0			,	0.0			,	3.0 / 10.0,
									0.0			,	6.0 / 10.0	,	0.0,
									0.0			,	0.0			,	6.0 / 10.0,
									0.0			,	10.0 / 10.0	,	10.0 / 10.0;
*/

	//ALL TRANSITIONS
	updateWeightMatrix <<	2.0 / 10.0	,	1.5 / 10.0	,	2.0 / 10.0,
												4.0 / 10.0	,	3.0 / 10.0	,	4.0 / 10.0,
												5.0 / 10.0	,	3.5 / 10.0	,	0.0 / 10.0,
												6.0 / 10.0	,	4.0 / 10.0	,	0.0 / 10.0,
												8.0 / 10.0	,	0.0			,	0.0,
												0.0			,	6.0 / 10.0	,	0.0, 
												10.0 / 10.0	,	0.0			,	0.0,
												0.0			,	0.0			,	6.0 / 10.0,
												0.0			,	8.0 / 10.0	,	0.0,
												0.0			,	0.0			,	8.0 / 10.0,
												0.0			,	10.0 / 10.0	,	10.0 / 10.0;

/*
	//ONLY Z<->W2<->W4
	updateWeightMatrix <<	2.5 / 10.0	,	2.0 / 10.0	,	2.5 / 10.0,
									5.0 / 10.0	,	4.0 / 10.0	,	5.0 / 10.0,
									6.5 / 10.0	,	5.0 / 10.0	,	6.5 / 10.0,
									8.0 / 10.0	,	6.0 / 10.0	,	8.0 / 10.0,
									10.0 / 10.0	,	0.0			,	0.0,
									0.0			,	8.0 / 10.0	,	0.0,
									0.0			,	0.0			,	0.0,
									0.0			,	0.0			,	0.0,
									0.0			,	10.0 / 10.0	,	0.0,
									0.0			,	0.0			,	10.0 / 10.0,
									0.0			,	0.0			,	0.0;
*/
/*
	//ONLY W2<->Z<->W4
	updateWeightMatrix <<	2.0 / 10.0	,	2.5 / 10.0	,	2.5 / 10.0,
									4.0 / 10.0	,	5.0 / 10.0	,	5.0 / 10.0,
									5.0 / 10.0	,	6.5 / 10.0	,	6.5 / 10.0,
									6.0 / 10.0	,	8.0 / 10.0	,	8.0 / 10.0,
									8.0 / 10.0	,	0.0			,	0.0,
									0.0			,	10.0 / 10.0	,	0.0,
									10.0 / 10.0	,	0.0			,	0.0,
									0.0			,	0.0			,	10.0 / 10.0,
									0.0			,	0.0			,	0.0,
									0.0			,	0.0			,	0.0,
									0.0			,	0.0			,	0.0;
*/
/*
	//ONLY Z<->W2
	updateWeightMatrix <<	2.5 / 10.0	,	2.0 / 10.0	,	0.0,
									5.0 / 10.0	,	4.0 / 10.0	,	0.0,
									6.5 / 10.0	,	5.0 / 10.0	,	0.0,
									8.0 / 10.0	,	6.0 / 10.0	,	0.0,
									10.0 / 10.0	,	0.0			,	0.0,
									0.0			,	8.0 / 10.0	,	0.0,
									0.0			,	0.0			,	0.0,
									0.0			,	0.0			,	0.0,
									0.0			,	0.0			,	0.0,
									0.0			,	0.0			,	0.0,
									0.0			,	10.0 / 10.0	,	0.0;
*/
/*
	//ONLY Z
	updateWeightMatrix <<	2.0 / 4.0,	0.0		,	0.0,
									4.0 / 4.0,	0.0		,	0.0,
									3.0 / 4.0,	0.0		,	0.0,
									4.0 / 4.0,	0.0		,	0.0,
									0.0		,	0.0		,	0.0,
									0.0		,	0.0		,	0.0, 
									0.0		,	0.0		,	0.0,
									0.0		,	0.0		,	0.0,
									0.0		,	0.0		,	0.0,
									0.0		,	0.0		,	0.0,
									0.0		,	0.0		,	0.0;
*/
	acceptedUpdates = matrix_t::Zero(nUpdateType, nStateType);
	proposedUpdates = matrix_t::Zero(nUpdateType, nStateType);
}

void mc::PrintAcceptanceMatrix(std::ostream& out)
{
	out << "Acceptance of updates:" << std::endl;
	for (uint_t i = 0; i < nUpdateType; ++i)
	{
		for (uint_t j = 0; j < nStateType; ++j)
		{
			if (proposedUpdates(i, j) == 0)
				out << 0 << " ";
			else
				out << acceptedUpdates(i, j) / proposedUpdates(i, j) << " ";
		}
		out << std::endl;
	}
}

void mc::do_update()
{
	if (!isInitialized)
	{
		std::cout << "Build G0 look up table";
		std::cout.flush();
		std::string filename = path + "g0lookup/rhom-L" + std::to_string(L);
		configSpace.BuildG0LookUpTable(filename);
		std::cout << "Done." << std::endl;
		configSpace.updateHandler.Init();
		if (!is_thermalized())
			std::cout << "Thermalization" << std::endl;
		isInitialized = true;
	}
	
	for (uint_t i = 0; i < nThermStep; ++i)
	{
		value_t r = rng();
		StateType state = configSpace.State();
		
		if (r < updateWeightMatrix(UpdateType::AddVertex, state))
		{
			bool result;
			if (state == StateType::Z)
				result = configSpace.AddRandomVertices<1, 0>();
			else if (state == StateType::W2)
				result = configSpace.AddRandomVertices<1, 1>();
			else if (state == StateType::W4)
				result = configSpace.AddRandomVertices<1, 2>();
			if (result)
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::AddVertex, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(UpdateType::AddVertex, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::RemoveVertex, state))
		{
			bool result;
			if (state == StateType::Z)
				result = configSpace.RemoveRandomVertices<1, 0>();
			else if (state == StateType::W2)
				result = configSpace.RemoveRandomVertices<1, 1>();
			else if (state == StateType::W4)
				result = configSpace.RemoveRandomVertices<1, 2>();
			if (result)
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::RemoveVertex, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(UpdateType::RemoveVertex, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::AddTwoVertices, state))
		{
			bool result;
			if (state == StateType::Z)
				result = configSpace.AddRandomVertices<2, 0>();
			else if (state == StateType::W2)
				result = configSpace.AddRandomVertices<2, 1>();
			else if (state == StateType::W4)
				result = configSpace.AddRandomVertices<2, 2>();
			if (result)
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::AddTwoVertices, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(UpdateType::AddTwoVertices, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::RemoveTwoVertices, state))
		{
			bool result;
			if (state == StateType::Z)
				result = configSpace.RemoveRandomVertices<2, 0>();
			else if (state == StateType::W2)
				result = configSpace.RemoveRandomVertices<2, 1>();
			else if (state == StateType::W4)
				result = configSpace.RemoveRandomVertices<2, 2>();
			if (result)
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::RemoveTwoVertices, state) += 1.0;
				++rebuildCnt;
			}
			proposedUpdates(UpdateType::RemoveTwoVertices, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::ZtoW2, state))
		{
			uint_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
			value_t preFactor = configSpace.lattice->Sites() * m * configSpace.beta * configSpace.zeta2;
			if (configSpace.AddRandomWorms<1, 0>(preFactor))
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::ZtoW2, state) += 1.0;
				configSpace.state = StateType::W2;
			}
			proposedUpdates(UpdateType::ZtoW2, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::W2toZ, state))
		{
			uint_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
			value_t preFactor = 1.0 / (configSpace.lattice->Sites() * m * configSpace.beta * configSpace.zeta2);
			if (configSpace.RemoveRandomWorms<1, 1>(preFactor))
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::W2toZ, state) += 1.0;
				configSpace.state = StateType::Z;
			}
			proposedUpdates(UpdateType::W2toZ, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::ZtoW4, state))
		{
			uint_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
			value_t preFactor = configSpace.lattice->Sites() * m * m * m * configSpace.beta * configSpace.zeta4;
			if (configSpace.AddRandomWorms<2, 0>(preFactor))
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::ZtoW4, state) += 1.0;
				configSpace.state = StateType::W4;
			}
			proposedUpdates(UpdateType::ZtoW4, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::W4toZ, state))
		{
			uint_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
			value_t preFactor = 1.0 / (configSpace.lattice->Sites() * m * m * m * configSpace.beta * configSpace.zeta4);
			if (configSpace.RemoveRandomWorms<2, 2>(preFactor))
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::W4toZ, state) += 1.0;
				configSpace.state = StateType::Z;
			}
			proposedUpdates(UpdateType::W4toZ, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::W2toW4, state))
		{
			uint_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
			value_t preFactor = (configSpace.lattice->Sites() * m * configSpace.zeta4) / configSpace.zeta2;
			if (configSpace.AddRandomWorms<1, 1>(preFactor))
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::W2toW4, state) += 1.0;
				configSpace.state = StateType::W4;
			}
			proposedUpdates(UpdateType::W2toW4, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::W4toW2, state))
		{
			uint_t m = configSpace.lattice->NeighborhoodCount(configSpace.nhoodDist);
			value_t preFactor = configSpace.zeta2 / (configSpace.lattice->Sites() * m * configSpace.zeta4);
			if (configSpace.RemoveRandomWorms<1, 2>(preFactor))
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::W4toW2, state) += 1.0;
				configSpace.state = StateType::W2;
			}
			proposedUpdates(UpdateType::W4toW2, state) += 1.0;
		}
		else if (r < updateWeightMatrix(UpdateType::shiftWorm, state))
		{
			bool result;
			if (state == StateType::Z)
				result = configSpace.ShiftWorm<0>();
			else if (state == StateType::W2)
				result = configSpace.ShiftWorm<1>();
			else if (state == StateType::W4)
				result = configSpace.ShiftWorm<2>();
			if (result)
			{
				configSpace.updateList.back() += " - success.";
				acceptedUpdates(UpdateType::shiftWorm, state) += 1.0;
			}
			proposedUpdates(UpdateType::shiftWorm, state) += 1.0;
		}

		value_t avgError = 0.0;
		value_t relError = 0.0;
			
		if (rebuildCnt == nRebuild)
		{
			//double cond = configSpace.updateHandler.StabilizeInvG(avgError, relError);
			double cond = configSpace.updateHandler.StabilizeInvG();
			measure.add("avgInvGError", avgError);
			measure.add("relInvGError", relError);
			measure.add("condition", cond);
			rebuildCnt = 0;
		}
	}
	
	if (!is_thermalized())
	{
		if ((sweep+1) % (nThermalize / 3) == 0)
		{
			std::cout << ".";
			std::cout.flush();
		}
	}
	if (sweep + 1 == nThermalize)
		std::cout << std::endl;
	++sweep;
}

void mc::do_measurement()
{	
	if ((sweep - nThermalize + 1) % (nMeasurements / 3) == 0)
	{
		std::cout << "Vertices: " << configSpace.updateHandler.GetVertexHandler().Vertices() << std::endl;
		std::cout << "Worms: " << configSpace.updateHandler.GetVertexHandler().Worms() << std::endl;
	}
	
	measure.add("<w>", configSpace.updateHandler.GetVertexHandler().Worms());
	uint_t R = 0;
	value_t sign, c;
	std::fill(corrVector.begin(), corrVector.end(), 0.0);
	switch (configSpace.State())
	{
		case StateType::Z:
			measure.add("deltaZ", 1.0);
			measure.add("deltaW2", 0.0);
			measure.add("deltaW4", 0.0);
			measure.add("k", configSpace.updateHandler.GetVertexHandler().Vertices());
			GetWithDef(exporderHistZ, configSpace.updateHandler.GetVertexHandler().Vertices(), 0) += 1;
			break;

		case StateType::W2:
			measure.add("deltaZ", 0.0);
			measure.add("deltaW2", 1.0);
			measure.add("deltaW4", 0.0);
			measure.add("k", 0.0);
			GetWithDef(exporderHistW2, configSpace.updateHandler.GetVertexHandler().Vertices(), 0) += 1;

			R = configSpace.updateHandler.GetVertexHandler().WormDistance();
			sign = configSpace.updateHandler.GetVertexHandler().WormParity();
			corrVector[R] = sign / configSpace.lattice->DistanceHistogram(R);
			break;
		case StateType::W4:
			measure.add("deltaZ", 0.0);
			measure.add("deltaW2", 0.0);
			measure.add("deltaW4", 1.0);
			measure.add("k", 0.0);
			break;
	}
	measure.add("Corr", corrVector);
}