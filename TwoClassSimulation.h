#pragma once 

#include "Simulation.h"

class TwoClassSimulation : public AbstractSimulation {
public:
	explicit TwoClassSimulation(const std::string& dataFilename, Toggles _toggles);
	double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) const override;
	OptionalPolicy heuristicPolicy() const override;
	OptionalPolicy optimalPolicy() const override;
	DistParameters eosConsumptionDistributionParametersForCustomer(int j, int u) const override;
};

