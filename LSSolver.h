#pragma once

#include <localsolver.h>
#include <vector>

class AbstractSimulation;

class RevenueComputationNativeFunction : public localsolver::LSNativeFunction {
public:
	RevenueComputationNativeFunction(AbstractSimulation &_sim, std::vector<std::vector<int>> &_scenarios) : sim(_sim), scenarios(_scenarios) {}
	virtual ~RevenueComputationNativeFunction() {}

	localsolver::lsdouble call(const localsolver::LSNativeContext& context) override;

private:
	AbstractSimulation &sim;
	std::vector<std::vector<int>> &scenarios;
};

void solveWithLS(AbstractSimulation &sim, std::vector<std::vector<int>> &scenarios);