#include <iostream>

#include "Evaluator.h"
#include "Simulation.h"
#include "Helpers.h"
#include "Stopwatch.h"

using namespace std;

Result AbstractEvaluator::computeOpt(const ResultList& results, bool printOpts) {
	Result optResult;
	optResult.profit = std::numeric_limits<float>::min();

	for(Result res : results) {
		if(res.profit >= optResult.profit) {
			optResult = res;
		}
	}

	if (printOpts) {
		cout << "Primary opt: " << optResult.toString() << endl;
		for (Result res : results) {
			if (res.profit == optResult.profit && res.bookingLimits != optResult.bookingLimits) {
				cout << "Alternative opt:" << res.toString() << endl;
			}
		}
	}

	return optResult;
}

Result AbstractEvaluator::solve(AbstractSimulation::ScenarioList& scenarios) {
	Stopwatch sw;
	sw.start();
	auto res = collectResults(scenarios);
	auto opt = computeOpt(res, false /*true*/);
	//cout << "Time passed = " << sw.lookAndReset() << endl;
	return opt;
}

ResultList Evaluator2D::collectResults(AbstractSimulation::ScenarioList &scenarios) {
	ResultList results(sim.getC() + 1);
	vector<int> bookingLimits(2);
	bookingLimits[0] = sim.getC();
	for(bookingLimits[1] = 0; bookingLimits[1] <= sim.getC(); bookingLimits[1]++) {
		results[bookingLimits[1]].bookingLimits = bookingLimits;
		results[bookingLimits[1]].profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));		
	}
	return results;
}

ResultList Evaluator3D::collectResults(AbstractSimulation::ScenarioList &scenarios) {
	ResultList results((int)(0.5 * (double)(sim.getC() + 1) * (double)(sim.getC() + 2)));
	vector<int> bookingLimits(sim.getNumClasses());
	bookingLimits[0] = sim.getC();
	int ctr = 0;
	for(bookingLimits[1] = 0; bookingLimits[1] <= bookingLimits[0]; bookingLimits[1]++) {
		for(bookingLimits[2] = 0; bookingLimits[2] <= bookingLimits[1]; bookingLimits[2]++) {
			results[ctr].bookingLimits = bookingLimits;
			results[ctr++].profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));			
		}
	}
	return results;
}
