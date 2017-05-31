//
// Created by André Schnabel on 23.05.17.
//

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <vector>
#include <iostream>
#include <array>

#include "Runner.h"
#include "Simulation.h"
#include "Evaluator.h"
#include "LSSolver.h"
#include "GurobiSolver.h"
#include "PSSolver.h"
#include "Helpers.h"

using namespace std;

namespace algo = boost::algorithm;
namespace fs = boost::filesystem;

struct Config {
	string instanceName, solverName;
	int numScenarios;

	Config(const string &instanceName, const string &solverName, int numScenarios) : instanceName(instanceName),
																					 solverName(solverName),
																					 numScenarios(numScenarios) {}
};

Config processArguments(const list<string> &args) {
	string	instanceName = "multi_data.json",
			solverName = "Gurobi",
			numScenarios = "150";

	auto assignFromArg = [](const string &prefix, const string &arg, string &out) {
		static auto getRhs = [](const string &line) {
			vector<string> parts;
			algo::split(parts, line, boost::is_any_of("="));
			assert(parts.size() == 2);
			return *next(parts.begin());
		};

		if(algo::starts_with(arg, prefix+"="))
			out = getRhs(arg);
	};

	for(string arg : args) {
		assignFromArg("instance", arg, instanceName);
		assignFromArg("solver", arg, solverName);
		assignFromArg("nscenarios", arg, numScenarios);
	}

	return { instanceName, solverName, stoi(numScenarios) };
}

map<string, function<BookingLimitOptimizer*()>> generateSolverNameToObjectMapping(const AbstractSimulation &sim) {
	return {
			{ "Gurobi", [&sim]() { return new GurobiOptimizer(sim); }},
			{ "LocalSolver", [&sim]() { return new LSOptimizer(sim); }},
			{ "ParticleSwarm", [&sim]() { return new PSSolver(sim); }},
			{ "FullEnumeration", [&sim]() { return new EvaluatorMultiDimensional(sim); }}
	};
};

void Runner::commandLine(const list<string> &args) {
	Config cfg = processArguments(args);
	MultiClassSimulation sim(cfg.instanceName+".json");

	auto scenarios = sim.generateScenarios(cfg.numScenarios, 42, AbstractSimulation::SamplingType::Descriptive);

	auto solverNameToObject = generateSolverNameToObjectMapping(sim);
	BookingLimitOptimizer *optimizer = solverNameToObject[cfg.solverName]();
	auto result = optimizer->solve(scenarios);
	cout << "Result = " << result.toString() << endl;
	delete optimizer;
}

list<string> instancesInDirectory(const string &dir) {
	list<string> lst;
	fs::path p(dir);
	fs::directory_iterator end;
	for(fs::directory_iterator it(p); it != end; ++it) {
		fs::path entry = it->path();
		if(fs::is_regular_file(entry) && algo::ends_with(entry.string(), ".json")) {
			string instanceName = entry.stem().string();
			lst.push_back(instanceName);
		}
	}
	return lst;
}

void Runner::benchmark(const string &dir) {
	static std::array<string, 1> solverNames = { "Gurobi" };
	//static std::array<string, 4> solverNames = { "Gurobi", "LocalSolver", "ParticleSwarm", "FullEnumeration" };

	static auto constructBookingLimitsCaption = [](int numClasses) {
		vector<string> bookingLimitStrs((unsigned long) numClasses);
		for(int i=0; i<numClasses; i++)
			bookingLimitStrs[i] = "b" + to_string((i+1));
		return algo::join(bookingLimitStrs, ";");
	};

	static auto constructBookingLimitsString = [](const vector<int> &bookingLimits) {
		vector<string> bookingLimitStrs(bookingLimits.size());
		for(int i=0; i<bookingLimits.size(); i++)
			bookingLimitStrs[i] = to_string(bookingLimits[i]);
		return algo::join(bookingLimitStrs, ";");
	};

	string bookingLimitsCaption = constructBookingLimitsCaption(3);

	for(const string &solverName : solverNames) {
		Helpers::spit("instance;profit;" + bookingLimitsCaption + "\n", solverName + "Results.txt");
	}

	list<string> instances = instancesInDirectory(dir);

	for(const string& instanceName : instances) {
		MultiClassSimulation sim(dir + "/" + instanceName + ".json");

		auto scenarios = sim.generateScenarios(150, 42, AbstractSimulation::SamplingType::Descriptive);// AbstractSimulation::SamplingType::Random);

		auto solverNameToObject = generateSolverNameToObjectMapping(sim);

		for(const string &solverName : solverNames) {
			auto optimizer = solverNameToObject[solverName]();
			auto result = optimizer->solve(scenarios);
			string bookingLimitsStr = constructBookingLimitsString(result.bookingLimits);
			Helpers::spitAppend(instanceName + ";" + to_string(result.profit) + ";" + bookingLimitsStr + "\n", solverName + "Results.txt");
			delete optimizer;
		}
	}
}
