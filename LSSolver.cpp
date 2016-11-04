#include "LSSolver.h"
#include "Simulation.h"
#include "Helpers.h"
#include "Evaluator.h"
#include <iostream>

using namespace std;
using namespace localsolver;

lsdouble RevenueComputationNativeFunction::call(const LSNativeContext& context) {
	int nclasses = sim.getNumClasses();

	vector<int> bookingLimits(nclasses);	
	for (int i = 0; i < nclasses; i++) {
		bookingLimits[i] = (int)context.getIntValue(i);
	}

	vector<double> profits(scenarios.size());
	for(int i=0; i<scenarios.size(); i++) {
		profits[i] = sim.objective(scenarios[i], bookingLimits);
	}

	return Helpers::vecAverage(profits);
}

void solveWithLS(AbstractSimulation& sim, AbstractSimulation::ScenarioList &scenarios) {
	RevenueComputationNativeFunction rfunc(sim, scenarios);

	LocalSolver ls;
	LSModel model = ls.getModel();

	vector<LSExpression> bookingLimits(sim.getNumClasses());

	auto obj = model.call(model.createNativeFunction(&rfunc));	

	for(int i=0; i<sim.getNumClasses(); i++) {
		bookingLimits[i] = model.intVar(0, sim.getC());
		obj.addOperand(bookingLimits[i]);
	}

	for (int i = 0; i<sim.getNumClasses()-1; i++) {
		model.constraint(bookingLimits[i] >= bookingLimits[i + 1]);
	}

	model.addObjective(obj, OD_Maximize);
	model.close();

	auto lsphase = ls.createPhase();
	lsphase.setTimeLimit(3);

	ls.solve();

	auto sol = ls.getSolution();

	AbstractEvaluator::Result res(sim.getNumClasses());
	res.profit = sol.getDoubleValue(obj);
	for (int i = 0; i < sim.getNumClasses(); i++)
		res.bookingLimits[i] = (int)sol.getIntValue(bookingLimits[i]);

	cout << "LocalSolver solution with native function: " << endl << res << endl;
}
