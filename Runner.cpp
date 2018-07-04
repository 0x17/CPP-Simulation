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
#include "EasyCSV.h"
#include "MultiClassSimulation.h"
#include "Helpers.h"
#include "Globals.h"

using namespace std;

namespace algo = boost::algorithm;
namespace fs = boost::filesystem;

double globals::timeLimit = 3.0;
bool globals::tracingEnabled = false;

struct CommandlineArguments {
	string instanceName, solverName;
	int numScenarios;
	bool tracingEnabled;
	double timeLimit;

	CommandlineArguments(
		const string &instanceName,
		const string &solverName,
		const int numScenarios,
		bool tracingEnabled,
		double timeLimit)
			: instanceName(instanceName),
			  solverName(solverName),
			  numScenarios(numScenarios),
			tracingEnabled(tracingEnabled),
			timeLimit(timeLimit) {}
};

void showUsage() {
	cout << "Usage: CPP-Simulation instance=multi_data solver=Gurobi nscenarios=150 timelimit=50.0 trace=true" << endl;
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

	CommandlineArguments cargs = { "multi_data.json", "Gurobi", 150, false, 30.0 };
	string numScenariosStr = to_string(cargs.numScenarios);
	string timeLimitStr = to_string(cargs.timeLimit);
	string tracingEnabledStr = to_string(cargs.tracingEnabled);

	for(const string &arg : args) {
		assignFromArg("instance", arg, cargs.instanceName);
		assignFromArg("solver", arg, cargs.solverName);
		assignFromArg("nscenarios", arg, numScenariosStr);
		assignFromArg("timelimit", arg, timeLimitStr);
		assignFromArg("trace", arg, tracingEnabledStr);
	}

	cargs.numScenarios = stoi(numScenariosStr);
	cargs.timeLimit = stod(timeLimitStr);
	cargs.tracingEnabled = tracingEnabledStr == "true";

	return cargs;
}

map<string, function<BookingLimitOptimizer*()>> generateSolverNameToObjectMapping(const AbstractSimulation &sim) {
	return {
			{ "Gurobi", [&sim]() { return new GurobiOptimizer(sim); }},
			{ "LocalSolver", [&sim]() { return new LSOptimizer(sim); }},
			{ "ParticleSwarm", [&sim]() { return new PSSolver(sim); }},
			{ "FullEnumeration", [&sim]() { return new EvaluatorMultiDimensional(sim); }}
	};
}



void Runner::commandLine(const list<string> &args) {
	const CommandlineArguments cfg = processArguments(args);

	globals::tracingEnabled = cfg.tracingEnabled;
	globals::timeLimit = cfg.timeLimit;

	const Toggles toggles("toggles.json");
	const MultiClassSimulation sim(cfg.instanceName+".json", toggles);
	const string actualInstanceName = fs::path(cfg.instanceName).stem().string();

	const auto scenarios = sim.generateDemandScenarios(cfg.numScenarios, 42, SamplingType::Descriptive);
	//const auto scenarios = sim.generateDemandScenariosBinomial(cfg.numScenarios, 42);
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

	const string ofn = cfg.solverName + "Results.txt";
	if (!fs::exists(ofn)) Helpers::spit("instanceName;profit\n", ofn);
	Helpers::spitAppend(actualInstanceName + ";" + to_string(result.profit) + "\n", ofn);
	//fs::rename(cfg.solverName + "Trace.txt", actualInstanceName + "_" + cfg.solverName + "Trace.txt");

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
	//static std::array<string, 1> solverNames = { "Gurobi" };
	//static std::array<string, 4> solverNames = { "Gurobi", "LocalSolver", "ParticleSwarm", "FullEnumeration" };
	static std::array<string, 3> solverNames = { "Gurobi", "LocalSolver", "ParticleSwarm" };

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

	MultiClassSimulation sim("multi_data_big.json", toggles);
	//MultiClassSimulation sim("multi_data.json", toggles);
	//MultiClassSimulation sim("vonfolie.json");
	//MultiClassSimulation sim("Instances/pinstance1.json");
    //MultiClassSimulation sim("multi_data_cvar_harder.json", toggles);

	sim.setRiskAversionParameters(0.8, 0.5);

	EvaluatorMultiDimensional evl(sim);
	LSOptimizer ls(sim);
	GurobiOptimizer gurobi(sim);
	PSSolver ps(sim);

	vector<BookingLimitOptimizer *> optimizers = { &evl, &ls, &gurobi, &ps };
	//vector<BookingLimitOptimizer *> optimizers = { &gurobi };

	const auto scenarios = sim.generateDemandScenarios(ntries, 42, SamplingType::Descriptive);

	list<pair<string, Result>> results;

	for (auto optimizer : optimizers) {
        cout << "Running: " << optimizer->getName() << endl;
		results.emplace_back(optimizer->getName(), optimizer->solve(scenarios, {}));
	}

	for (auto res : results) {
		cout << endl << res.first << " results:" << endl;
		cout << res.second.toString() << endl << endl;
		cout << "Comparison objective: " << sim.averageRevenueOfSimulation(res.second.bookingLimits, scenarios) << endl;
	}
}

void serializeDemandScenariosToDisk(const DemandScenarioList &scenarios, const std::string &outFilename) {
	const auto headerParts = Helpers::constructVector<string>(scenarios.getN(), [](int j) { return "j" + to_string(j+1);  });
    EasyCSV ecsv("scenario;" + boost::algorithm::join(headerParts, ";"));

	for (int s = 0; s < scenarios.getM(); s++) {
		const auto scenario = scenarios.row(s);
		ecsv.addRow(Helpers::constructVector<string>(scenarios.getN() + 1, [s, &scenario](int j) { return to_string(j == 0 ? s : scenario[j - 1]); }));
	}
	ecsv.persist(outFilename);
}

Matrix<Result> computeAndPersistAlphaPsiVariations(AbstractSimulation &sim, const DemandScenarioList &scenarios, BookingLimitOptimizer &blo) {
	auto persistScenarioProfits = [&scenarios, &sim](int alphaCtr, int psiCtr, const vector<int> &bookingLimits) {
		EasyCSV profitScenarios("scenario;profit");
		for (int sindex = 0; sindex < scenarios.getM(); sindex++) {
			const auto demands = scenarios.row(sindex);
			profitScenarios.addRow(vector<double>{ (double)sindex, sim.objective(demands, bookingLimits) });
		}
		profitScenarios.persist("psi" + to_string(psiCtr) + "_alpha" + to_string(alphaCtr) + "_profits");
	};

	Matrix<Result> resmx(11, 11, {});
	EasyCSV ecsv("alpha;psi;profit;b1;b2;b3");

	int psiCtr = 0;
	for (double psi = 0.0; psi <= 1.0; psi += 0.1, psiCtr++) {
		int alphaCtr = 0;
		for (double alpha = 0.0; alpha <= 1.0; alpha += 0.1, alphaCtr++) {
			sim.setRiskAversionParameters(alpha, psi);
			auto res = blo.solve(scenarios, {});
			resmx(alphaCtr, psiCtr) = res;
			ecsv.addRow({alpha, psi, res.profit, (double) res.bookingLimits[0], (double) res.bookingLimits[1], (double) res.bookingLimits[2]});
			persistScenarioProfits(alphaCtr, psiCtr, res.bookingLimits);
		}
	}

	ecsv.persist("alphaPsiVariations");
	return resmx;
}

void persistAlphaPsiVariationsTable(Matrix<Result> &resmx) {
	const auto colCaptions = Helpers::constructVector<std::string>(11, [](int i) { return to_string(i * 0.1); });
	EasyCSV ecsv("TABLE;"+ boost::algorithm::join(colCaptions, ";"));

	auto resultCellStr = [](Result &res) {
		const vector<double> parts = {res.profit, (double) res.bookingLimits[0], (double) res.bookingLimits[1],
									  (double) res.bookingLimits[2]};
		const auto sparts = Helpers::constructVector<string>(parts.size(),
															 [&parts](int i) { return to_string(parts[i]); });
		return boost::algorithm::join(sparts, ",");
		//return to_string(res.profit);
	};

	for (int row = 0; row < 11; row++) {
		const vector<string> parts = Helpers::constructVector<string>(12, [&resultCellStr, &resmx, row](int col) {
			return (col == 0) ? to_string(row * 0.1) : resultCellStr(resmx(row, col - 1));
		});
		ecsv.addRow(parts);
	}

	ecsv.persist("alphaPsiVariationsTable");
}

void Runner::alphaPsiVariations() {
	const Toggles toggles;
	MultiClassSimulation sim("multi_data_cvar.json", toggles);
	const unique_ptr<BookingLimitOptimizer> blo = make_unique<GurobiOptimizer>(sim);

    assert(sim.getNumClasses() == 3);

	const auto scenarios = sim.generateDemandScenarios(150, 42, SamplingType::Descriptive);
    serializeDemandScenariosToDisk(scenarios, "scenariosDescriptive");

	serializeDemandScenariosToDisk(sim.generateDemandScenarios(150, 42, SamplingType::Random), "scenarios");

	Matrix<Result> resmx = computeAndPersistAlphaPsiVariations(sim, scenarios, *blo);
	persistAlphaPsiVariationsTable(resmx);
}
