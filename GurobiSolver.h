#pragma once

#include "Simulation.h"

class GurobiOptimizer : public BookingLimitOptimizer {
public:
	explicit GurobiOptimizer(const AbstractSimulation& _sim) : BookingLimitOptimizer("Gurobi", _sim) {}
	Result solve(const DemandScenarioList& scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) override;

private:
	template<class Func>
	Result solveCommon(const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc, Func modelBuilder);
};