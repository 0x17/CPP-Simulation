#include "LSSolver.h"
#include "Simulation.h"
#include "Helpers.h"

using namespace std;
using namespace localsolver;

lsdouble RevenueComputationNativeFunction::call(const LSNativeContext& context) {
	static Helpers::Tracer tr("LocalSolverTrace");
	static double bks = numeric_limits<double>::lowest();

	int nclasses = sim.getNumClasses();

	vector<int> bookingLimits(nclasses);	
	for (int i = 0; i < nclasses; i++) {
		bookingLimits[i] = (int)context.getIntValue(i);
	}

	double obj = sim.averageRevenueOfSimulation(bookingLimits, scenarios);
	if(obj > bks) {
		bks = obj;
	}

	tr.intervalTrace(bks);

	return obj;
}

LSOptimizer::LSOptimizer(AbstractSimulation& _sim): BookingLimitOptimizer("LocalSolverNative", _sim), rfunc(_sim), bookingLimits(sim.getNumClasses()) {
	LSModel model = ls.getModel();

	obj = model.call(model.createNativeFunction(&rfunc));

	for (int i = 0; i<sim.getNumClasses(); i++) {
		bookingLimits[i] = model.intVar(0, sim.getC());
		obj.addOperand(bookingLimits[i]);
	}

	model.constraint(bookingLimits[0] == sim.getC());
	for (int i = 0; i<sim.getNumClasses() - 1; i++) {
		model.constraint(bookingLimits[i] >= bookingLimits[i + 1]);
	}

	model.addObjective(obj, OD_Maximize);
	model.close();

	if (heuristicBookingLimits) {
		for (int j = 0; j<sim.getNumClasses(); j++) {
			bookingLimits[j].setIntValue((*heuristicBookingLimits)[j]);
		}
	}

	ls.getParam().setNbThreads(1);
}

Result LSOptimizer::solve(vector<vector<int>>& scenarios) {
	rfunc.setScenarios(scenarios);

	auto lsphase = ls.createPhase();
	lsphase.setTimeLimit(5);

	ls.solve();

	auto sol = ls.getSolution();

	Result res(sim.getNumClasses());
	res.profit = sol.getDoubleValue(obj);
	for (int i = 0; i < sim.getNumClasses(); i++)
		res.bookingLimits[i] = (int)sol.getIntValue(bookingLimits[i]);

	return res;
}
