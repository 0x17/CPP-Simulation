#pragma once

#include <vector>
#include "Simulation.h"

class GurobiOptimizer : public BookingLimitOptimizer {
public:
	explicit GurobiOptimizer(const AbstractSimulation& _sim) : BookingLimitOptimizer("Gurobi", _sim) {}
	Result solve(const ScenarioList& scenarios) override;

private:
	Result solveWithEconomiesOfScale(const ScenarioList &scenarios);
	Result solveWithNewFormulation(const ScenarioList &scenarios);
	Result solveWithOldFormulation(const ScenarioList &scenarios);
};