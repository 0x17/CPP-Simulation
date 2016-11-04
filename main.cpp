#include <iostream>

#include "Simulation.h"
#include "Evaluator.h"
#include "LSSolver.h"
#include "GurobiSolver.h"

using namespace std;

const int ntries = 50;

void runOptimizers() {
	cout << "Number of scenarios: " << ntries << endl << endl;

	MultiClassSimulation sim;

	Evaluator3D evl(sim);
	LSOptimizer ls(sim);
	GurobiOptimizer gurobi(sim);

	vector<BookingLimitOptimizer *> optimizers = {  &evl, &ls, &gurobi };

	auto scenarios = sim.generateScenarios(ntries);

	for(auto optimizer : optimizers) {
		auto res = optimizer->solve(scenarios);
		cout << endl << optimizer->getName() << " results:" << endl;
		cout << res << endl << endl;
	}
}

int main() {
	runOptimizers();	
	getchar();
    return 0;
}