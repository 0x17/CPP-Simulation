//
// Created by André Schnabel on 30.10.16.
//

#include <algorithm>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include "json11.hpp"

#include "Simulation.h"
#include "Helpers.h"
#include "Globals.h"

using namespace std;

string Result::toString() const {
	bool first = true;
    string out = "profit=" + to_string(profit) + ",bookingLimits=(";
	for (double bl : this->bookingLimits) {
		if (!first) out += ", ";
		out += to_string(bl);
		first = false;
	}
	out += ")";
    return out;
}

AbstractSimulation::AbstractSimulation(const string &dataFilename) {
    string errMsg;
    json11::Json obj = json11::Json::parse(Helpers::slurp(dataFilename), errMsg);
	if(!errMsg.empty()) {
		throw runtime_error("Unable to parse " + dataFilename + " error = " + errMsg);
	}
    C = obj["capacity"].int_value();
	for(const auto &c : obj["clients"].array_items()) {
		customers.emplace_back(c);
	}
    numClasses = (int)customers.size();
}

vector<double> AbstractSimulation::runSimulation(const vector<int> &bookingLimits, const ScenarioList &scenarios) const {
    vector<double> revenues(static_cast<unsigned long>(scenarios.getM()));
    for(int i=0; i<revenues.size(); i++) {
        revenues[i] = objective(scenarios.row(i), bookingLimits);
    }
    return revenues;
}

double AbstractSimulation::averageRevenueOfSimulation(const std::vector<int>& bookingLimits, const ScenarioList& scenarios) const {
	return Helpers::vecAverage(runSimulation(bookingLimits, scenarios));
}

vector<double> AbstractSimulation::statisticalMeansOfScenarios(ScenarioList& scenarios) {
	vector<double> customerMeans(static_cast<unsigned long>(scenarios.getN()), 0.0);
	int scenarioCount = scenarios.getM();

	scenarios.foreach([&customerMeans](int i, int j, int val) {
		customerMeans[j] += val;
	});

	for (double &customerMean : customerMeans) {
		customerMean /= (double)scenarioCount;
	}

	return customerMeans;
}

vector<double> AbstractSimulation::statisticalStandardDeviationsOfScenarios(ScenarioList& scenarios) {
	auto means = statisticalMeansOfScenarios(scenarios);
	vector<double> stddevs(static_cast<unsigned long>(scenarios.getN()), 0.0);
	int scenarioCount = scenarios.getM();

	scenarios.foreach([&stddevs, &means](int i, int j, int val) {
		stddevs[j] += pow(val - means[j], 2);
	});

	for (double &stddev : stddevs) {
		stddev = sqrt(stddev / (double)scenarioCount);
	}

	return stddevs;
}

double TwoClassSimulation::objective(const vector<int>& demands, const vector<int>& bookingLimits) const {
	int n2 = (int)floor(min((double)bookingLimits[1], demands[1] * customers[1].consumptionPerReq) / customers[1].consumptionPerReq);
	int n1 = (int)floor(min(demands[0] * customers[0].consumptionPerReq, C - n2 * customers[1].consumptionPerReq));
	return n1 * customers[0].revenuePerReq + n2 * customers[1].revenuePerReq;
}

OptionalPolicy TwoClassSimulation::heuristicPolicy() const {
	vector<int> bookingLimits(2);
	bookingLimits[0] = (int)floor(min(customers[0].expD * customers[0].consumptionPerReq, (double)C) / customers[0].consumptionPerReq);
	bookingLimits[1] = (int)floor(min(customers[1].expD * customers[1].consumptionPerReq, (C-bookingLimits[0]*customers[0].consumptionPerReq)) / customers[1].consumptionPerReq);
	return bookingLimits;
}

OptionalPolicy TwoClassSimulation::optimalPolicy() const {	
	if (customers[0].consumptionPerReq == 1 && customers[1].consumptionPerReq == 1) {
		double x = (customers[0].revenuePerReq - customers[1].revenuePerReq) / customers[0].revenuePerReq;
		vector<int> bookingLimits = { C, C - (int)floor(Helpers::invNormal(x, customers[0].expD, customers[0].devD)) };
		return bookingLimits;
	}
	return OptionalPolicy();
}

double MultiClassSimulation::objective(const vector<int>& demands, const vector<int>& bookingLimits) const {
	double	residualCapacity = C,
			profit = 0.0;

	for(int ix = numClasses - 1; ix >= 0; ix--) {
		int n;
		if (globals::ECONOMY_OF_SCALE_ENABLED) {
			double blCapacityLeft = (bookingLimits[ix] - (C - residualCapacity));
			int fittingInBL = 0;
			for (int k = 1; k <= ceil(blCapacityLeft); k++) {
				if (k * eosConsumption(ix, k) <= blCapacityLeft)
					fittingInBL = k;
			}
			n = min(demands[ix], fittingInBL);
			residualCapacity -= n * eosConsumption(ix, n);
			//cout << "n" << ix << "=" << n << endl;
		}
		else {
			n = min(demands[ix], (int)floor((bookingLimits[ix] - (C - residualCapacity)) / customers[ix].consumptionPerReq));
			residualCapacity -= n * customers[ix].consumptionPerReq;
		}

		profit += n * customers[ix].revenuePerReq;
	}

	//cout << "Profit == " << profit << endl;

	return profit;
}

OptionalPolicy MultiClassSimulation::heuristicPolicy() const {
	vector<int> bookingLimits((unsigned long)numClasses);
	bookingLimits[0] = C;	
	for(int j=1; j<numClasses; j++) {
		bookingLimits[j] = (int)floor( min(customers[j].expD * customers[j].consumptionPerReq, (double)bookingLimits[j-1]) );
	}
	return bookingLimits;
}

OptionalPolicy MultiClassSimulation::optimalPolicy() const {
	return OptionalPolicy();
}

double MultiClassSimulation::eosConsumption(int j, int u) const {
	double deltaX = customers[j].expD;
	double deltaY = customers[j].consumptionPerReq * 0.2;
	return max(1.0, customers[j].consumptionPerReq - (boost::math::erf(4 * u / deltaX - 2) * deltaY / 2.0 + deltaY / 2.0));
}

Scenario AbstractSimulation::pickDemands(int scenarioIx, int numScenarios) {
	Scenario demands((unsigned long)numClasses);
	int customerIx = 0;
	for(Customer &c : customers) {
		double samplingResult = Helpers::pickNormal(c.expD, c.devD);
		demands[customerIx++] = max(0, (int)round(samplingResult));
	}
	return demands;
}

Scenario AbstractSimulation::pickDemandsDescriptive(int scenarioIx, int numScenarios, LUTList &lutList) {
    Scenario demands((unsigned long)numClasses);
    int customerIx = 0;
    for(Customer &c : customers) {
		double samplingResult = Helpers::pickNextWithLUT(lutList[customerIx]);
		demands[customerIx++] = max(0, (int)round(samplingResult));
    }
    return demands;
}

ScenarioList AbstractSimulation::generateScenarios(int ntries, int seed, SamplingType stype) {
	Helpers::resetSeed(seed);

	ScenarioList scenarios(ntries, static_cast<int>(customers.size()));

	if(stype == SamplingType::Descriptive) {
		LUTList lookupTables(customers.size());
		for (int i = 0; i < customers.size(); i++) {
			lookupTables[i] = Helpers::generateNormalDistributionDescriptiveSamplingLUT(ntries, customers[i].expD, customers[i].devD);
		}
		for(int i=0; i<ntries; i++) {
			scenarios.setRow(i, pickDemandsDescriptive(i, ntries, lookupTables));
		}
	} else {
		for(int i=0; i<ntries; i++) {
			scenarios.setRow(i, pickDemands(i, ntries));
		}
	}

    return scenarios;
}

Customer::Customer(const json11::Json &obj) :
        name(obj["name"].string_value()),
        expD(obj["expD"].number_value()),
        devD(obj["devD"].number_value()),
        description(obj["description"].string_value()),
        consumptionPerReq(obj["consumptionPerReq"].number_value()),
        revenuePerReq(obj["revenuePerReq"].number_value())
{}
