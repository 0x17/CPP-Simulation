#pragma once

#include <localsolver.h>
#include <vector>
#include "Simulation.h"

class RevenueComputationNativeFunction : public localsolver::LSNativeFunction {
public:
	RevenueComputationNativeFunction(const AbstractSimulation &_sim) : sim(_sim), scenarios(emptyScenarios) {}
	RevenueComputationNativeFunction(const AbstractSimulation &_sim, DemandScenarioList &_scenarios) : sim(_sim), scenarios(_scenarios) {}
	virtual ~RevenueComputationNativeFunction() {}

	localsolver::lsdouble call(const localsolver::LSNativeContext& context) override;

	void setScenarios(const DemandScenarioList &_scenarios) { this->scenarios = _scenarios; }

private:
	const AbstractSimulation &sim;
	DemandScenarioList &scenarios;
	DemandScenarioList emptyScenarios = {};
};

class LSOptimizer : public BookingLimitOptimizer {
public:
	LSOptimizer(const AbstractSimulation& _sim);
	Result solve(const DemandScenarioList& scenarios, const boost::optional<ConsumptionScenarioList> &consumptionScenarios) override;
private:
	localsolver::LocalSolver ls;
	RevenueComputationNativeFunction rfunc;
	std::vector<localsolver::LSExpression> bookingLimits;
	localsolver::LSExpression obj;
};