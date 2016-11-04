#include <iostream>

#include "Simulation.h"
#include "Evaluator.h"
#include "LSSolver.h"
#include "GurobiSolver.h"
#include "Helpers.h"

using namespace std;

void runOptimizers() {
	const int ntries = 10;

	cout << "Number of scenarios: " << ntries << endl << endl;

	MultiClassSimulation sim;

	Evaluator3D evl(sim);
	LSOptimizer ls(sim);
	GurobiOptimizer gurobi(sim);

	vector<BookingLimitOptimizer *> optimizers = {  &evl, &ls, &gurobi };

	auto scenarios = sim.generateScenarios(ntries, AbstractSimulation::SamplingType::Descriptive);

	for(auto optimizer : optimizers) {
		auto res = optimizer->solve(scenarios);
		cout << endl << optimizer->getName() << " results:" << endl;
		cout << res << endl << endl;
	}
}

void effectOfDescriptiveSampling() {
	MultiClassSimulation sim;
	Evaluator3D evl(sim);

	auto refScenarios = sim.generateScenarios(10000000, AbstractSimulation::SamplingType::Random);

	Helpers::spit("ntries;profitRandom;profitDescriptive\n", "profitfornscen.txt");

	for(int ntries = 1; ntries <= 50; ntries += 1) {
		auto scenariosRand = sim.generateScenarios(ntries, AbstractSimulation::SamplingType::Random);
		auto resRand = evl.solve(scenariosRand);
		auto profitRandomSampling = Helpers::vecAverage(sim.runSimulation(resRand.bookingLimits, refScenarios));

		auto scenariosDescr = sim.generateScenarios(ntries, AbstractSimulation::SamplingType::Descriptive);
		auto resDescr = evl.solve(scenariosDescr);
		auto profitDescriptiveSampling = Helpers::vecAverage(sim.runSimulation(resDescr.bookingLimits, refScenarios));

		Helpers::spitAppend(to_string(ntries) + ";" + to_string(profitRandomSampling) + ";" + to_string(profitDescriptiveSampling) + "\n", "profitfornscen.txt");
	}
}

int main() {
	//runOptimizers();
	effectOfDescriptiveSampling();
	getchar();
    return 0;
}