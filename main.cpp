#include <iostream>
#include <boost/algorithm/string.hpp>

#include "Simulation.h"
#include "Evaluator.h"
#include "LSSolver.h"
#include "GurobiSolver.h"
#include "PSSolver.h"
#include "Helpers.h"
#include "Runner.h"
#include <thread>

using namespace std;

void runOptimizers();

int main(int argc, const char **argv) {
	//runOptimizers();
	//testInverseNormal();
	//effectOfDescriptiveSampling();

	/*cout << "Press [Return] to quit!" << endl;
	getchar();*/

	Runner::commandLine(Helpers::extractArguments(argc, argv));

	return 0;
}

void runOptimizers() {
	const int ntries = 150;

	cout << "Number of scenarios: " << ntries << endl << endl;

	//MultiClassSimulation sim("multi_data_big.json");
	MultiClassSimulation sim("multi_data.json");

	EvaluatorMultiDimensional evl(sim);
	LSOptimizer ls(sim);
	GurobiOptimizer gurobi(sim);
	PSSolver ps(sim);

	vector<BookingLimitOptimizer *> optimizers = {  &evl, &ls, &gurobi, &ps };
	//vector<BookingLimitOptimizer *> optimizers = { &ls };

	auto scenarios = sim.generateScenarios(ntries, 42, AbstractSimulation::SamplingType::Descriptive);

	list<pair<string, Result>> results;

	for(auto optimizer : optimizers) {
		results.push_back(make_pair(optimizer->getName(), optimizer->solve(scenarios)));
	}

	for(auto res : results) {
		cout << endl << res.first << " results:" << endl;
		cout << res.second.toString() << endl << endl;
	}
}

