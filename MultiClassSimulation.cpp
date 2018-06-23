
#include <boost/math/special_functions.hpp>
#include "MultiClassSimulation.h"

using namespace std;


MultiClassSimulation::MultiClassSimulation(const std::string& dataFilename, Toggles _toggles) : AbstractSimulation(dataFilename, _toggles) {
}

double MultiClassSimulation::objective(const vector<int>& demands, const vector<int>& bookingLimits) const {
	double	cumulativeConsumption = 0,
		profit = 0.0;

	for (int ix = numClasses - 1; ix >= 0; ix--) {
		int n;
		if (toggles.economyOfScale) {
			const double blCapacityLeft = (bookingLimits[ix] - cumulativeConsumption);
			int fittingInBL = 0;
			for (int k = 1; k <= ceil(blCapacityLeft); k++) {
				if (k * eosConsumption(ix, k) <= blCapacityLeft)
					fittingInBL = k;
			}
			n = min(demands[ix], fittingInBL);
			cumulativeConsumption += n * eosConsumption(ix, n);
		}
		else {
			n = min(demands[ix], (int)floor((bookingLimits[ix] - cumulativeConsumption) / customers[ix].consumptionPerReqMean));
			cumulativeConsumption += n * customers[ix].consumptionPerReqMean;
		}

		profit += n * customers[ix].revenuePerReq;
	}

	return profit;
}

OptionalPolicy MultiClassSimulation::heuristicPolicy() const {
	vector<int> bookingLimits((unsigned long)numClasses);
	bookingLimits[0] = C;
	for (int j = 1; j<numClasses; j++) {
		bookingLimits[j] = (int)floor(min(customers[j].expD * customers[j].consumptionPerReqMean, (double)bookingLimits[j - 1]));
	}
	return bookingLimits;
}

OptionalPolicy MultiClassSimulation::optimalPolicy() const {
	return OptionalPolicy();
}

double MultiClassSimulation::eosConsumption(int j, int u) const {
	const double deltaX = customers[j].expD;
	const double deltaY = customers[j].consumptionPerReqMean * 0.2;
	return max(1.0, customers[j].consumptionPerReqMean - (boost::math::erf(4 * u / deltaX - 2) * deltaY / 2.0 + deltaY / 2.0));
}

DistParameters MultiClassSimulation::eosConsumptionDistributionParametersForCustomer(int j, int u) const {
	const double mu = eosConsumption(j, u);
	const double sigma = mu / 10.0;
	return DistParameters{ mu, sigma };
}