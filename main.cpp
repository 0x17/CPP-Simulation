#include <iostream>

#include "Simulation.h"
#include "Evaluator.h"
#include "LSSolver.h"
#include "GurobiSolver.h"
#include "Helpers.h"

#include <boost/algorithm/string.hpp>
#include "PSSolver.h"

using namespace std;

void runOptimizers() {
	const int ntries = 10;

	cout << "Number of scenarios: " << ntries << endl << endl;

	MultiClassSimulation sim;

	Evaluator3D evl(sim);
	LSOptimizer ls(sim);
	GurobiOptimizer gurobi(sim);
	PSSolver ps(sim);

	vector<BookingLimitOptimizer *> optimizers = {  &evl, &ls, &gurobi, &ps };

	auto scenarios = sim.generateScenarios(ntries, 42, AbstractSimulation::SamplingType::Descriptive);

	for(auto optimizer : optimizers) {
		auto res = optimizer->solve(scenarios);
		cout << endl << optimizer->getName() << " results:" << endl;
		cout << res.toString() << endl << endl;
	}
}

inline string dot2comma(string s) {
	return boost::algorithm::replace_all_copy(s, ".", ",");
}

void effectOfDescriptiveSampling() {
	MultiClassSimulation sim;
	Evaluator3D evl(sim);

	auto refScenarios = sim.generateScenarios(100000, 23, AbstractSimulation::SamplingType::Descriptive);

	//auto means = AbstractSimulation::statisticalMeansOfScenarios(refScenarios);
	//auto stddevs = AbstractSimulation::statisticalStandardDeviationsOfScenarios(refScenarios);

	Helpers::spit("ntries;profitRandom;profitDescriptive\n", "profitfornscen.txt");
	Helpers::spit("blimitsRand;blimitsDescriptive\n", "blimits.txt");
	Helpers::spit("diffMeanRand;diffMeanDescr;diffStdDevRand;diffStdDevDescr\n", "diffmeans.txt");

	for(int ntries = 1; ntries <= 2000; ntries += 50) {
		cout << "n=" << ntries << endl;
		auto scenariosRand = sim.generateScenarios(ntries, 42, AbstractSimulation::SamplingType::Random);
		auto resRand = evl.solve(scenariosRand);
		auto profitRandomSampling = Helpers::vecAverage(sim.runSimulation(resRand.bookingLimits, refScenarios));

		auto scenariosDescr = sim.generateScenarios(ntries, 42, AbstractSimulation::SamplingType::Descriptive);
		auto resDescr = evl.solve(scenariosDescr);
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

int main() {
	runOptimizers();
	//effectOfDescriptiveSampling();
	//getchar();
    return 0;
}