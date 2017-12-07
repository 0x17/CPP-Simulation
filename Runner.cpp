//
// Created by André Schnabel on 23.05.17.
//

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <iostream>

#include "Runner.h"
#include "Simulation.h"
#include "Evaluator.h"
#include "LSSolver.h"
#include "GurobiSolver.h"
#include "PSSolver.h"
#include "Globals.h"

using namespace std;

namespace algo = boost::algorithm;
namespace fs = boost::filesystem;

struct CommandlineArguments {
	string instanceName, solverName;
	int numScenarios;

	CommandlineArguments(const string &instanceName,
		   const string &solverName,
		   const int numScenarios)
			: instanceName(instanceName),
			  solverName(solverName),
			  numScenarios(numScenarios) {}
};

void showUsage() {
	cout << "Usage: CPP-Simulation instance=multi_data.json solverName=Gurobi numScenarios=150" << endl;
	cout << "Solver names = { Gurobi, LocalSolver, ParticleSwarm, FullEnumeration }" << endl;
}

string getRighthandSide(const string &line) {
	vector<string> parts;
	algo::split(parts, line, boost::is_any_of("="));
	assert(parts.size() == 2);
	return *next(parts.begin());
};

void assignFromArg(const string &prefix, const string &arg, string &out) {
	if (algo::starts_with(arg, prefix + "="))
		out = getRighthandSide(arg);
};

CommandlineArguments processArguments(const list<string> &args) {
	if(args.empty()) {
		showUsage();
		throw runtime_error("Not enough args!");
	}

	CommandlineArguments cargs = { "multi_data.json", "Gurobi", 150 };
	string numScenariosStr = to_string(cargs.numScenarios);

	for(const string &arg : args) {
		assignFromArg("instance", arg, cargs.instanceName);
		assignFromArg("solver", arg, cargs.solverName);
		assignFromArg("nscenarios", arg, numScenariosStr);
	}

	cargs.numScenarios = stoi(numScenariosStr);
	return cargs;
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
	CommandlineArguments cfg = processArguments(args);
	const Toggles toggles("toggles.json");
	const MultiClassSimulation sim(cfg.instanceName+".json", toggles);

	const auto scenarios = sim.generateDemandScenarios(cfg.numScenarios, 42, SamplingType::Descriptive);
	auto solverNameToObject = generateSolverNameToObjectMapping(sim);
	BookingLimitOptimizer *optimizer = solverNameToObject[cfg.solverName]();

	Result result;
	if(toggles.stochasticConsumptions) {
		ConsumptionScenarioFunc csfunc = sim.generateConsumptionScenarioFunc(cfg.numScenarios, 42, SamplingType::Descriptive, toggles.economyOfScale ? -1 : 0);
		const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc = csfunc;
		result = optimizer->solve(scenarios, consumptionScenarioFunc);
	} else {
		result = optimizer->solve(scenarios, {});
	}
	cout << "Result = " << result.toString() << endl;
	delete optimizer;
}

list<string> instancesInDirectory(const string &dir) {
	list<string> lst;
	fs::path p(dir);
	const fs::directory_iterator end;
	for(fs::directory_iterator it(p); it != end; ++it) {
		fs::path entry = it->path();
		if(fs::is_regular_file(entry) && algo::ends_with(entry.string(), ".json")) {
			const string instanceName = entry.stem().string();
			lst.push_back(instanceName);
		}
	}
	return lst;
}

void Runner::benchmark(const string &dir) {
	const Toggles toggles("Toggles.json");
	static std::array<string, 1> solverNames = { "Gurobi" };
	//static std::array<string, 4> solverNames = { "Gurobi", "LocalSolver", "ParticleSwarm", "FullEnumeration" };

	static auto constructBookingLimitsCaption = [](int numClasses) {
		vector<string> bookingLimitStrs(static_cast<unsigned long>(numClasses));
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

	const string bookingLimitsCaption = constructBookingLimitsCaption(3);

	for(const string &solverName : solverNames) {
		Helpers::spit("instance;profit;" + bookingLimitsCaption + "\n", solverName + "Results.txt");
	}

	list<string> instances = instancesInDirectory(dir);

	for(const string& instanceName : instances) {
		MultiClassSimulation sim(dir + "/" + instanceName + ".json", toggles);

		const auto scenarios = sim.generateDemandScenarios(150, 42, SamplingType::Descriptive);// AbstractSimulation::SamplingType::Random);

		auto solverNameToObject = generateSolverNameToObjectMapping(sim);

		for(const string &solverName : solverNames) {
			auto optimizer = solverNameToObject[solverName]();
			auto result = optimizer->solve(scenarios, {});
			const string bookingLimitsStr = constructBookingLimitsString(result.bookingLimits);
			Helpers::spitAppend(instanceName + ";" + to_string(result.profit) + ";" + bookingLimitsStr + "\n", solverName + "Results.txt");
			delete optimizer;
		}
	}
}

void Runner::runOptimizers() {
	const Toggles toggles("Toggles.json");
	const int ntries = 150;

	cout << "Number of scenarios: " << ntries << endl << endl;

	//MultiClassSimulation sim("multi_data_big.json");
	MultiClassSimulation sim("multi_data.json", toggles);
	//MultiClassSimulation sim("vonfolie.json");
	//MultiClassSimulation sim("Instances/pinstance1.json");

	EvaluatorMultiDimensional evl(sim);
	LSOptimizer ls(sim);
	GurobiOptimizer gurobi(sim);
	PSSolver ps(sim);

	vector<BookingLimitOptimizer *> optimizers = { &evl, &ls, &gurobi, &ps };
	//vector<BookingLimitOptimizer *> optimizers = { &gurobi };

	const auto scenarios = sim.generateDemandScenarios(ntries, 42, SamplingType::Descriptive);

	list<pair<string, Result>> results;

	for (auto optimizer : optimizers) {
		results.emplace_back(optimizer->getName(), optimizer->solve(scenarios, {}));
	}

	for (auto res : results) {
		cout << endl << res.first << " results:" << endl;
		cout << res.second.toString() << endl << endl;
		cout << "Comparison objective: " << sim.averageRevenueOfSimulation(res.second.bookingLimits, scenarios) << endl;
	}
}

void Runner::alphaPsiVariations() {
	const Toggles toggles;
	MultiClassSimulation sim("multi_data.json", toggles);
	unique_ptr<BookingLimitOptimizer> blo = make_unique<GurobiOptimizer>(sim);

	const auto scenarios = sim.generateDemandScenarios(150, 42, SamplingType::Descriptive);

	const string OUT_FN = "alphaPsiVariations.csv";
	Helpers::spit("alpha;psi;profit;b1;b2;b3\n", OUT_FN);

	int alphaCtr = 0;
	int psiCtr = 0;	

	Matrix<Result> resmx(11, 11, {});
	
	for(double alpha = 0.0; alpha <= 1.0; alpha += 0.1, alphaCtr++) {
		for(double psi = 0.0; psi <= 1.0; psi += 0.1, psiCtr++) {
			
			sim.setRiskAversionParameters(alpha, psi);
			auto res = blo.get()->solve(scenarios, {});

			resmx(alphaCtr, psiCtr) = res;

			const vector<double> parts = { alpha, psi, res.profit, (double)res.bookingLimits[0], (double)res.bookingLimits[1], (double)res.bookingLimits[2] };
			const auto sparts = Helpers::constructVector<string>(parts.size(), [&parts](int i) { return to_string(parts[i]); });
			const string resline = boost::algorithm::join(sparts, ";");
			Helpers::spitAppend(resline, OUT_FN);
		}
	}

	const string OUT_FN2 = "alphaPsiVariationsTable.csv";
	const auto colCaptions = Helpers::constructVector<std::string>(11, [](int i) { return to_string(i*0.1); });
	const string colNamesStr = boost::algorithm::join(colCaptions, ";");
	Helpers::spit("TABLE;"+colNamesStr+"\n", OUT_FN2);

	auto resultCellStr = [](Result &res) {
		const vector<double> parts = { res.profit, (double)res.bookingLimits[0], (double)res.bookingLimits[1], (double)res.bookingLimits[2] };
		const auto sparts = Helpers::constructVector<string>(parts.size(), [&parts](int i) { return to_string(parts[i]); });
		return boost::algorithm::join(sparts, ",");
	};

	for(int row = 0; row < 11; row++) {
		const vector<string> parts = Helpers::constructVector<string>(12, [&resultCellStr, &resmx, row](int col) {
			return (col == 0) ? to_string(col*0.1) : resultCellStr(resmx(row, col - 1));
		});
		const string rowline = boost::algorithm::join(parts, ";");
		Helpers::spitAppend(rowline, OUT_FN2);
	}
}
