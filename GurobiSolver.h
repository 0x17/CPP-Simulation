#pragma once

#include <vector>
#include "Simulation.h"

class GurobiOptimizer : public BookingLimitOptimizer {
public:
	explicit GurobiOptimizer(const AbstractSimulation& _sim) : BookingLimitOptimizer("Gurobi", _sim) {}
	Result solve(std::vector<std::vector<int>>& scenarios) override;
};