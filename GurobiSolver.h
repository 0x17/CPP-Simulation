#pragma once

#include <vector>
#include "Simulation.h"

class GurobiOptimizer : public BookingLimitOptimizer {
public:
	explicit GurobiOptimizer(const AbstractSimulation& _sim) : BookingLimitOptimizer("Gurobi", _sim) {}
	Result solve(const DemandScenarioList& scenarios) override;

private:
	template<class Func>
	Result solveCommon(const DemandScenarioList &scenarios, Func modelBuilder);
};