#pragma once

#include <vector>
#include "Simulation.h"

class GurobiOptimizer : public BookingLimitOptimizer {
public:
	explicit GurobiOptimizer(const AbstractSimulation& _sim) : BookingLimitOptimizer("Gurobi", _sim) {}
	Result solve(std::vector<std::vector<int>>& scenarios) override;

private:
	Result solveWithEconomiesOfScale(std::vector<std::vector<int>>& scenarios);
	Result solveWithNewFormulation(std::vector<std::vector<int>>& scenarios);
	Result solveWithOldFormulation(std::vector<std::vector<int>>& scenarios);
};