#include <iostream>

#include "Evaluator.h"
#include "Simulation.h"
#include "Helpers.h"
#include "Stopwatch.h"
#include <boost/algorithm/string/join.hpp>

using namespace std;

Result AbstractEvaluator::computeOpt(const ResultList& results, bool printOpts) {
	Result optResult;
	optResult.profit = numeric_limits<float>::lowest();

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

ResultList Evaluator2D::collectResults(AbstractSimulation::ScenarioList &scenarios) const {
	ResultList results(sim.getC() + 1);
	vector<int> bookingLimits(2);
	bookingLimits[0] = sim.getC();
	for(bookingLimits[1] = 0; bookingLimits[1] <= sim.getC(); bookingLimits[1]++) {
		results[bookingLimits[1]].bookingLimits = bookingLimits;
		results[bookingLimits[1]].profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));		
	}
	return results;
}

ResultList EvaluatorMultiDimensional::collectResults(AbstractSimulation::ScenarioList &scenarios) const {
	ResultList resultList(1);
	vector<int> bookingLimits(sim.getNumClasses());
	Helpers::Tracer tr("FullEnumerationTrace");
	Stopwatch sw;
	const int timelimit = 30;

	resultList[0].profit = numeric_limits<double>::lowest();

	bookingLimits[0] = sim.getC();

	sw.start();
	double tstart = sw.look();

	function<void(int)> recursiveCollector = [&](int classIndex) {
		if (timelimit != -1 && sw.look() - tstart >= (double)timelimit * 1000.0)
			return;

		if (classIndex == bookingLimits.size()) {
			tr.intervalTrace(resultList[0].profit);

			double obj = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));
			if(obj > resultList[0].profit) {
				resultList[0].profit = obj;
				resultList[0].bookingLimits = bookingLimits;
			}
		}
		else {
			for (bookingLimits[classIndex] = 0; bookingLimits[classIndex] <= bookingLimits[classIndex - 1]; bookingLimits[classIndex]++) {
				recursiveCollector(classIndex + 1);
			}
		}
	};

	recursiveCollector(1);

	return resultList;
}

ResultList Evaluator3D::collectResults(AbstractSimulation::ScenarioList &scenarios) const {
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
