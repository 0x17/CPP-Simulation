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
	if(errMsg.size() > 0) {
		throw new runtime_error("Unable to parse " + dataFilename + " error = " + errMsg);
	}
    C = obj["capacity"].int_value();
	for(auto c : obj["clients"].array_items()) {
		customers.push_back(Customer(c));
	}
    numClasses = (int)customers.size();
}

vector<double> AbstractSimulation::runSimulation(const vector<int> &bookingLimits, ScenarioList &scenarios) const {
    vector<double> revenues(scenarios.size());
    for(int i=0; i<scenarios.size(); i++) {
        revenues[i] = objective(scenarios[i], bookingLimits);
    }
    return revenues;
}

double AbstractSimulation::averageRevenueOfSimulation(const std::vector<int>& bookingLimits, ScenarioList& scenarios) const {
	return Helpers::vecAverage(runSimulation(bookingLimits, scenarios));
}

vector<double> AbstractSimulation::statisticalMeansOfScenarios(ScenarioList& scenarios) {
	vector<double> means = Helpers::constructVector<double>((int)scenarios[0].size(), [&](int ix) { return 0.0; });
	for(auto scenario : scenarios) {
		for(int i=0; i<means.size(); i++) {
			means[i] += scenario[i];
		}
	}
	for(int i=0; i<means.size(); i++) {
		means[i] /= scenarios.size();
	}
	return means;
}

vector<double> AbstractSimulation::statisticalStandardDeviationsOfScenarios(ScenarioList& scenarios) {
	auto means = statisticalMeansOfScenarios(scenarios);
	vector<double> stddevs = Helpers::constructVector<double>((int)scenarios[0].size(), [&](int ix) { return 0.0; });
	for (auto scenario : scenarios) {
		for (int i = 0; i<stddevs.size(); i++) {
			stddevs[i] += pow(scenario[i] - means[i], 2);
		}
	}
	for (int i = 0; i<stddevs.size(); i++) {
		stddevs[i] = sqrt(stddevs[i] / scenarios.size());
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
		int n = min(demands[ix], (int)floor((bookingLimits[ix] - (C - residualCapacity)) / customers[ix].consumptionPerReq));
		residualCapacity -= n * eosConsumption(ix, n);
		profit += n * customers[ix].revenuePerReq;
	}

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
	if(globals::ECONOMY_OF_SCALE_ENABLED) {
		double deltaX = customers[j].expD;
		double deltaY = customers[j].consumptionPerReq * 0.2;
		return max(1.0, customers[j].consumptionPerReq - (boost::math::erf(4 * u / deltaX - 2) * deltaY / 2.0 + deltaY / 2.0));
	} else {
		return customers[j].consumptionPerReq;
	}
}

AbstractSimulation::Scenario AbstractSimulation::pickDemands(int scenarioIx, int numScenarios) {
	Scenario demands((unsigned long)numClasses);
	int customerIx = 0;
	for(Customer &c : customers) {
		double samplingResult = Helpers::pickNormal(c.expD, c.devD);
		demands[customerIx++] = max(0, (int)round(samplingResult));
	}
	return demands;
}

AbstractSimulation::Scenario AbstractSimulation::pickDemandsDescriptive(int scenarioIx, int numScenarios, LUTList &lutList) {
    Scenario demands((unsigned long)numClasses);
    int customerIx = 0;
    for(Customer &c : customers) {
		double samplingResult = Helpers::pickNextWithLUT(lutList[customerIx]);
		demands[customerIx++] = max(0, (int)round(samplingResult));
    }
    return demands;
}

AbstractSimulation::ScenarioList AbstractSimulation::generateScenarios(int ntries, int seed, SamplingType stype) {
	Helpers::resetSeed(seed);

	ScenarioList scenarios((unsigned long) ntries);

	if(stype == SamplingType::Descriptive) {
		LUTList lookupTables(customers.size());
		for (int i = 0; i < customers.size(); i++) {
			lookupTables[i] = Helpers::generateNormalDistributionDescriptiveSamplingLUT(ntries, customers[i].expD, customers[i].devD);
		}
		for(int i=0; i<ntries; i++) {
			scenarios[i] = pickDemandsDescriptive(i, ntries, lookupTables);
		}
	} else {
		for(int i=0; i<ntries; i++) {
			scenarios[i] = pickDemands(i, ntries);
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
