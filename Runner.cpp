//
// Created by André Schnabel on 23.05.17.
//

#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>
#include <iostream>

#include "Runner.h"
#include "Simulation.h"
#include "Evaluator.h"
#include "LSSolver.h"
#include "GurobiSolver.h"
#include "PSSolver.h"

using namespace std;
namespace algo = boost::algorithm;

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

void Runner::commandLine(const list<string> &args) {
	Config cfg = processArguments(args);
	MultiClassSimulation sim(cfg.instanceName);

	auto scenarios = sim.generateScenarios(cfg.numScenarios, 42, AbstractSimulation::SamplingType::Descriptive);

	map<string, function<BookingLimitOptimizer*()>> solverNameToObject = {
			{ "Gurobi", [&sim]() { return new GurobiOptimizer(sim); }},
			{ "LocalSolver", [&sim]() { return new LSOptimizer(sim); }},
			{ "ParticleSwarm", [&sim]() { return new PSSolver(sim); }},
			{ "FullEnumeration", [&sim]() { return new EvaluatorMultiDimensional(sim); }}
	};

	BookingLimitOptimizer *optimizer = solverNameToObject[cfg.solverName]();
	auto result = optimizer->solve(scenarios);
	cout << "Result = " << result.toString() << endl;
}