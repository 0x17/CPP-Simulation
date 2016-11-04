#pragma once

#include <localsolver.h>
#include <vector>
#include "Simulation.h"


class RevenueComputationNativeFunction : public localsolver::LSNativeFunction {
public:
	RevenueComputationNativeFunction(AbstractSimulation &_sim) : sim(_sim), scenarios(emptyScenarios) {}
	RevenueComputationNativeFunction(AbstractSimulation &_sim, std::vector<std::vector<int>> &_scenarios) : sim(_sim), scenarios(_scenarios) {}
	virtual ~RevenueComputationNativeFunction() {}

	localsolver::lsdouble call(const localsolver::LSNativeContext& context) override;

	void setScenarios(std::vector<std::vector<int>> &_scenarios) { this->scenarios = _scenarios; }

private:
	AbstractSimulation &sim;
	std::vector<std::vector<int>> &scenarios, emptyScenarios = {};
};

class LSOptimizer : public BookingLimitOptimizer {
public:
	LSOptimizer(AbstractSimulation& _sim);
	Result solve(std::vector<std::vector<int>>& scenarios) override;
private:
	localsolver::LocalSolver ls;
	RevenueComputationNativeFunction rfunc;
	std::vector<localsolver::LSExpression> bookingLimits;
	localsolver::LSExpression obj;
};