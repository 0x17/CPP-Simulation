#include <iostream>
#include <boost/algorithm/string/join.hpp>

#include "Evaluator.h"
#include "Helpers.h"
#include "Globals.h"

using namespace std;

Result AbstractEvaluator::extractOptimumFromList(const ResultList& results, bool printOpts) {
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

Result AbstractEvaluator::solve(const DemandScenarioList& scenarios) {
	Stopwatch sw;
	sw.start();
	/*auto res = collectResults(scenarios);
	auto opt = extractOptimumFromList(res, false);*/
	auto opt = computeOptimum(scenarios);
	cout << endl << "Time passed = " << sw.lookAndReset() << endl;
	return opt;
}

ResultList Evaluator2D::collectResults(const DemandScenarioList &scenarios) const {
	ResultList results(static_cast<unsigned long>(sim.getC() + 1));
	vector<int> bookingLimits(2);
	bookingLimits[0] = sim.getC();
	for(bookingLimits[1] = 0; bookingLimits[1] <= sim.getC(); bookingLimits[1]++) {
		results[bookingLimits[1]].bookingLimits = bookingLimits;
		results[bookingLimits[1]].profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));
		cout << "Evaluating booking limit b_1=" << bookingLimits[1] << "\r" << flush;
	}
	return results;
}

Result Evaluator2D::computeOptimum(const DemandScenarioList& scenarios) const {
	Result opt;
	vector<int> bookingLimits(2);
	bookingLimits[0] = sim.getC();
	for (bookingLimits[1] = 0; bookingLimits[1] <= sim.getC(); bookingLimits[1]++) {
		double profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));
		if(profit > opt.profit) {
			opt.bookingLimits = bookingLimits;
			opt.profit = profit;
		}
		cout << "Evaluating booking limit b_1=" << bookingLimits[1] << "\r" << flush;
	}
	return opt;
}

ResultList Evaluator3D::collectResults(const DemandScenarioList &scenarios) const {
	ResultList results(static_cast<unsigned long>((int)(0.5 * (double)(sim.getC() + 1) * (double)(sim.getC() + 2))));
	vector<int> bookingLimits(static_cast<unsigned long>(sim.getNumClasses()));
	bookingLimits[0] = sim.getC();
	int ctr = 0;
	for (bookingLimits[1] = 0; bookingLimits[1] <= bookingLimits[0]; bookingLimits[1]++) {
		for (bookingLimits[2] = 0; bookingLimits[2] <= bookingLimits[1]; bookingLimits[2]++) {
			results[ctr].bookingLimits = bookingLimits;
			results[ctr++].profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));
			cout << "Evaluating booking limit b_1=" << bookingLimits[1] << ",b_2=" << bookingLimits[2] << "\r" << flush;
		}
	}
	return results;
}

Result Evaluator3D::computeOptimum(const DemandScenarioList& scenarios) const {
	Result opt;
	vector<int> bookingLimits(static_cast<unsigned long>(sim.getNumClasses()));
	bookingLimits[0] = sim.getC();
	for (bookingLimits[1] = 0; bookingLimits[1] <= bookingLimits[0]; bookingLimits[1]++) {
		for (bookingLimits[2] = 0; bookingLimits[2] <= bookingLimits[1]; bookingLimits[2]++) {
			double profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));;
			if(profit > opt.profit) {
				opt.bookingLimits = bookingLimits;
				opt.profit = profit;
			}			
			cout << "Evaluating booking limit b_1=" << bookingLimits[1] << ",b_2=" << bookingLimits[2] << "\r" << flush;
		}
	}
	return opt;
}

ResultList EvaluatorMultiDimensional::collectResults(const DemandScenarioList &scenarios) const {
	ResultList resultList(1);
	vector<int> bookingLimits(static_cast<unsigned long>(sim.getNumClasses()));
	Helpers::Tracer tr("FullEnumerationTrace");
	Stopwatch sw;
	const int timelimit = 15;

	resultList[0].profit = numeric_limits<double>::lowest();

	bookingLimits[0] = sim.getC();

	sw.start();
	double tstart = sw.look();

	function<void(int)> recursiveCollector = [&](int classIndex) {
		if (timelimit != -1 && sw.look() - tstart >= (double)timelimit * 1000.0)
			return;

		if (classIndex == bookingLimits.size()) {
			tr.intervalTrace(resultList[0].profit);

			double obj = sim.objectiveWithGlobalSettings(bookingLimits, scenarios);
			if(obj > resultList[0].profit) {
				resultList[0].profit = obj;
				resultList[0].bookingLimits = bookingLimits;
			}
		}
		else {
			for (bookingLimits[classIndex] = 0; bookingLimits[classIndex] <= bookingLimits[classIndex - 1]; bookingLimits[classIndex]++) {
				cout << "Booking limit b_" << classIndex << "=" << bookingLimits[classIndex] << "\r" << flush;
				recursiveCollector(classIndex + 1);
			}
		}
	};

	recursiveCollector(1);

	return resultList;
}

Result EvaluatorMultiDimensional::computeOptimum(const DemandScenarioList& scenarios) const {
	Result opt;
	vector<int> bookingLimits(static_cast<unsigned long>(sim.getNumClasses()));
	Helpers::Tracer tr("FullEnumerationTrace");
	Stopwatch sw;
	const double timelimit = globals::TIME_LIMIT;

	bookingLimits[0] = sim.getC();

	sw.start();
	double tstart = sw.look();

	function<void(int)> recursiveCollector = [&](int classIndex) {
		if (timelimit != -1.0 && sw.look() - tstart >= timelimit * 1000.0)
			return;

		if (classIndex == bookingLimits.size()) {
			tr.intervalTrace(opt.profit);

			double obj = sim.objectiveWithGlobalSettings(bookingLimits, scenarios);
			if (obj > opt.profit) {
				opt.profit = obj;
				opt.bookingLimits = bookingLimits;
			}
		}
		else {
			for (bookingLimits[classIndex] = 0; bookingLimits[classIndex] <= bookingLimits[classIndex - 1]; bookingLimits[classIndex]++) {
				cout << "Booking limit b_" << classIndex << "=" << bookingLimits[classIndex] << "\r" << flush;
				recursiveCollector(classIndex + 1);
			}
		}
	};

	recursiveCollector(1);

	return opt;
}
