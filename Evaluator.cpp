#include <iostream>

#include "Evaluator.h"
#include "Simulation.h"
#include "Helpers.h"
#include "Stopwatch.h"
#include <boost/algorithm/string/join.hpp>

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

ResultList EvaluatorMultiDimensional::collectResults(AbstractSimulation::ScenarioList &scenarios) {
	ResultList results((int)(0.5 * (double)(sim.getC() + 1) * (double)(sim.getC() + 2)));
	vector<int> bookingLimits(sim.getNumClasses());
	
	int ctr = 0;

	auto buildFirstRow = [&]() {
		std::vector<std::string> classCaptions(sim.getNumClasses() - 1);
		for (int i = 1; i < sim.getNumClasses(); i++) {
			classCaptions[i] = "b" + to_string(i + 1);
		}
		string firstRow = boost::algorithm::join(classCaptions, ";");
		firstRow += ";obj\n";
		return firstRow;
	};

	auto buildRow = [&]() {
		std::vector<std::string> blstrs(sim.getNumClasses() - 1);
		for (int i = 1; i<sim.getNumClasses(); i++) {
			blstrs[i] = to_string(bookingLimits[i]);
		}
		string row = boost::algorithm::join(blstrs, ";");
		row += ";" + to_string(results[ctr - 1].profit) + "\n";
		return row;
	};

	Helpers::spit(buildFirstRow(), "fullEnumObjectives.txt");

	bookingLimits[0] = sim.getC();

	std::function<void(int)> recursiveCollector = [&](int classIndex) {
		if (classIndex == bookingLimits.size()) {
			results[ctr].bookingLimits = bookingLimits;
			results[ctr++].profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));			
			Helpers::spitAppend(buildRow(), "fullEnumObjectives.txt");
		}
		else {
			for (bookingLimits[classIndex] = 0; bookingLimits[classIndex] <= bookingLimits[classIndex - 1]; bookingLimits[classIndex]++) {
				recursiveCollector(classIndex + 1);
			}
		}
	};

	recursiveCollector(1);

	/*for(bookingLimits[1] = 0; bookingLimits[1] <= bookingLimits[0]; bookingLimits[1]++) {
		for(bookingLimits[2] = 0; bookingLimits[2] <= bookingLimits[1]; bookingLimits[2]++) {
			results[ctr].bookingLimits = bookingLimits;
			results[ctr++].profit = Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));
			Helpers::spitAppend(to_string(bookingLimits[1]) + ";" + to_string(bookingLimits[2]) + ";" + to_string(results[ctr - 1].profit) + "\n", "fullEnumObjectives.txt");
		}
	}*/

	return results;
}
