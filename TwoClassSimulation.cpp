
#include <cmath>
#include "TwoClassSimulation.h"
#include "EasyCSV.h"
#include "Helpers.h"

using namespace std;


TwoClassSimulation::TwoClassSimulation(const std::string& dataFilename, Toggles _toggles) : AbstractSimulation(dataFilename, _toggles) {
}

double TwoClassSimulation::objective(const vector<int>& demands, const vector<int>& bookingLimits) const {
	int n2 = (int)floor(min((double)bookingLimits[1], demands[1] * customers[1].consumptionPerReqMean) / customers[1].consumptionPerReqMean);
	int n1 = (int)floor(min(demands[0] * customers[0].consumptionPerReqMean, C - n2 * customers[1].consumptionPerReqMean));
	return n1 * customers[0].revenuePerReq + n2 * customers[1].revenuePerReq;
}

OptionalPolicy TwoClassSimulation::heuristicPolicy() const {
	vector<int> bookingLimits(2);
	bookingLimits[0] = (int)floor(min(customers[0].expD * customers[0].consumptionPerReqMean, (double)C) / customers[0].consumptionPerReqMean);
	bookingLimits[1] = (int)floor(min(customers[1].expD * customers[1].consumptionPerReqMean, (C - bookingLimits[0] * customers[0].consumptionPerReqMean)) / customers[1].consumptionPerReqMean);
	return bookingLimits;
}

OptionalPolicy TwoClassSimulation::optimalPolicy() const {
	if (customers[0].consumptionPerReqMean == 1 && customers[1].consumptionPerReqMean == 1) {
		double x = (customers[0].revenuePerReq - customers[1].revenuePerReq) / customers[0].revenuePerReq;
		vector<int> bookingLimits = { C, C - (int)floor(Helpers::invNormal(x, customers[0].expD, customers[0].devD)) };
		return bookingLimits;
	}
	return OptionalPolicy();
}

DistParameters TwoClassSimulation::eosConsumptionDistributionParametersForCustomer(int u, int j) const {
	return DistParameters{ customers[j].consumptionPerReqMean, customers[j].consumptionPerReqStdDev };
}
