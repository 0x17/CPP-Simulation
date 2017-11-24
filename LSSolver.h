#pragma once

#include <localsolver.h>
#include <vector>
#include "Simulation.h"

class RevenueComputationNativeFunction : public localsolver::LSNativeFunction {
public:
	RevenueComputationNativeFunction(const AbstractSimulation &_sim) : sim(_sim), scenarios(emptyScenarios) {}
	RevenueComputationNativeFunction(const AbstractSimulation &_sim, ScenarioList &_scenarios) : sim(_sim), scenarios(_scenarios) {}
	virtual ~RevenueComputationNativeFunction() {}

	localsolver::lsdouble call(const localsolver::LSNativeContext& context) override;

	void setScenarios(const ScenarioList &_scenarios) { this->scenarios = _scenarios; }

private:
	const AbstractSimulation &sim;
	ScenarioList &scenarios;
	ScenarioList emptyScenarios = {};
};

class LSOptimizer : public BookingLimitOptimizer {
public:
	LSOptimizer(const AbstractSimulation& _sim);
	Result solve(const ScenarioList& scenarios) override;
private:
	localsolver::LocalSolver ls;
	RevenueComputationNativeFunction rfunc;
	std::vector<localsolver::LSExpression> bookingLimits;
	localsolver::LSExpression obj;
};