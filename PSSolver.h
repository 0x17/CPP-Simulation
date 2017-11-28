#pragma once

#include "Simulation.h"

class PSSolver : public BookingLimitOptimizer {
public:	
	explicit PSSolver(const AbstractSimulation &_sim);
	Result solve(const DemandScenarioList& scenarios, const boost::optional<ConsumptionScenarioList> &consumptionScenarios) override;
};

