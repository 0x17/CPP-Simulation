#pragma once

#include "Simulation.h"

class PSSolver : public BookingLimitOptimizer {
public:	
	explicit PSSolver(const AbstractSimulation &_sim);
	Result solve(const ScenarioList& scenarios) override;
};

