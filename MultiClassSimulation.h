#pragma once 

#include "Simulation.h"

class MultiClassSimulation : public AbstractSimulation {
public:
	explicit MultiClassSimulation(const std::string& dataFilename, Toggles _toggles);
	double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) const override;
	OptionalPolicy heuristicPolicy() const override;
	OptionalPolicy optimalPolicy() const override;
	double eosConsumption(int j, int u) const;
	DistParameters eosConsumptionDistributionParametersForCustomer(int j, int u) const override;
};

