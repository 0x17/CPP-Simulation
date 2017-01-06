#pragma once

#include "Simulation.h"

class PSSolver : public BookingLimitOptimizer {
public:	
	PSSolver(AbstractSimulation &_sim);
	Result solve(std::vector<std::vector<int>>& scenarios) override;
};

