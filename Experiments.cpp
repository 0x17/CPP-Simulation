//
// Created by André Schnabel on 23.05.17.
//

#include "Experiments.h"
#include "Simulation.h"
#include "Evaluator.h"
#include "Helpers.h"
#include "GurobiSolver.h"
#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <cmath>
using namespace std;

enum class SolutionMethod {
	FullEnumeration = 0,
	Gurobi
};

void Experiments::effectOfDescriptiveSampling() {
	auto dot2comma = [](string s) {
		return boost::algorithm::replace_all_copy(s, ".", ",");
	};

	//MultiClassSimulation sim("multi_data.json");
	//MultiClassSimulation sim("multi_data_alternative.json");
	MultiClassSimulation sim("multi_data_big.json");

	EvaluatorMultiDimensional evl(sim);

	auto refScenarios = sim.generateDemandScenarios(100000, 23, SamplingType::Descriptive);

	//auto means = AbstractSimulation::statisticalMeansOfScenarios(refScenarios);
	//auto stddevs = AbstractSimulation::statisticalStandardDeviationsOfScenarios(refScenarios);

	Helpers::spit("ntries;profitRandom;profitDescriptive\n", "profitfornscen.txt");
	Helpers::spit("blimitsRand;blimitsDescriptive\n", "blimits.txt");
	Helpers::spit("diffMeanRand;diffMeanDescr;diffStdDevRand;diffStdDevDescr\n", "diffmeans.txt");
	Helpers::spit("ntries;solvetime\n", "solvetimeforntries.txt");

	SolutionMethod method = SolutionMethod::Gurobi;

	for(int ntries = 1; ntries <= /*200*/ 106; ntries += 5) {
		cout << "n=" << ntries << endl;

		auto scenariosRand = sim.generateDemandScenarios(ntries, 42, SamplingType::Random);
		auto scenariosDescr = sim.generateDemandScenarios(ntries, 42, SamplingType::Descriptive);

		Result resRand, resDescr;
		switch(method) {
		case SolutionMethod::FullEnumeration:
			resRand = evl.solve(scenariosRand);
			resDescr = evl.solve(scenariosDescr);
			break;
		case SolutionMethod::Gurobi:
			GurobiOptimizer optimizer(sim);
			resRand = optimizer.solve(scenariosRand);
			resDescr = optimizer.solve(scenariosDescr);
			break;
		}

		auto profitRandomSampling = Helpers::vecAverage(sim.runSimulation(resRand.bookingLimits, refScenarios));
		auto profitDescriptiveSampling = Helpers::vecAverage(sim.runSimulation(resDescr.bookingLimits, refScenarios));

		Helpers::spitAppend(dot2comma(to_string(ntries) + ";" +to_string(profitRandomSampling) + ";" + to_string(profitDescriptiveSampling) + "\n"), "profitfornscen.txt");
		Helpers::spitAppend(resRand.toString() + ";" + resDescr.toString() + "\n", "blimits.txt");

		auto diffMeansRand = AbstractSimulation::statisticalMeansOfScenarios(scenariosRand);
		auto diffStddevsRand = AbstractSimulation::statisticalStandardDeviationsOfScenarios(scenariosRand);
		auto diffMeansDescr = AbstractSimulation::statisticalMeansOfScenarios(scenariosDescr);
		auto diffStddevsDescr = AbstractSimulation::statisticalStandardDeviationsOfScenarios(scenariosDescr);

		for (int i = 0; i < diffMeansRand.size(); i++) {
			diffMeansRand[i] = abs(diffMeansRand[i] - sim.getCustomer(i).expD);
			diffStddevsRand[i] = abs(diffStddevsRand[i] - sim.getCustomer(i).devD);
			diffMeansDescr[i] = abs(diffMeansDescr[i] - sim.getCustomer(i).expD);
			diffStddevsDescr[i] = abs(diffStddevsDescr[i] - sim.getCustomer(i).devD);
		}

		Helpers::spitAppend(dot2comma(
			to_string(Helpers::vecAverage(diffMeansRand))
			+ ";" + to_string(Helpers::vecAverage(diffMeansDescr))
			+ ";" + to_string(Helpers::vecAverage(diffStddevsRand))
			+ ";" + to_string(Helpers::vecAverage(diffStddevsDescr))
			+ "\n"), "diffmeans.txt");
	}
}

void Experiments::testInverseNormal() {
	const double epsilon = 0.00001;
	double inv = Helpers::invNormal(0.5, 0.0, 1.0);
	assert(abs(inv - 0.0) < epsilon);
}