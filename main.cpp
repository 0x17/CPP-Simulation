#include <iostream>

#include "Simulation.h"
#include "Evaluator.h"
#include "LSSolver.h"
#include "Stopwatch.h"
#include "GurobiSolver.h"

using namespace std;

int ntries = 10000;

void runSimulation3D() {
	cout << "Number of scenarios: " << ntries << endl << endl;
	cout << "Full enumeration results: " << endl;

	MultiClassSimulation sim;
	Evaluator3D evl(sim);
	Stopwatch sw;

	auto scenarios = sim.generateScenarios(ntries);

	sw.start();
	auto res = evl.collectResults(scenarios);
	auto opt = evl.computeOpt(res, true);
	cout << "Time passed = " << sw.lookAndReset() << endl;
	//cout << "Optimal solution from full enumeration: " << endl << opt << endl << endl;
	cout << endl << endl;

	solveWithLS(sim, scenarios);
	solveWithGurobi(sim, scenarios);
}

int main() {
	runSimulation3D();	
	getchar();
    return 0;
}