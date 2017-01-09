//
// Created by André Schnabel on 30.10.16.
//

#include <algorithm>
#include <cmath>
#include "json11.hpp"

#include "Simulation.h"
#include "Helpers.h"

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

vector<double> AbstractSimulation::runSimulation(const vector<int> &bookingLimits, ScenarioList &scenarios) {
    vector<double> revenues(scenarios.size());
    for(int i=0; i<scenarios.size(); i++) {
        revenues[i] = objective(scenarios[i], bookingLimits);
    }
    return revenues;
}

double AbstractSimulation::averageRevenueOfSimulation(const std::vector<int>& bookingLimits, ScenarioList& scenarios) {
	return Helpers::vecAverage(runSimulation(bookingLimits, scenarios));
}

std::vector<double> AbstractSimulation::statisticalMeansOfScenarios(ScenarioList& scenarios) {
	vector<double> means = Helpers::constructVector<double>(scenarios[0].size(), [&](int ix) { return 0.0; });
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

std::vector<double> AbstractSimulation::statisticalStandardDeviationsOfScenarios(ScenarioList& scenarios) {
	auto means = statisticalMeansOfScenarios(scenarios);
	vector<double> stddevs = Helpers::constructVector<double>(scenarios[0].size(), [&](int ix) { return 0.0; });
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

double TwoClassSimulation::objective(const vector<int>& demands, const vector<int>& bookingLimits) {
	int n2 = (int)floor(min(bookingLimits[1], demands[1] * customers[1].consumptionPerReq) / customers[1].consumptionPerReq);
	int n1 = (int)floor(min(demands[0] * customers[0].consumptionPerReq, C - n2 * customers[1].consumptionPerReq));
	return n1 * customers[0].revenuePerReq + n2 * customers[1].revenuePerReq;
}

OptionalPolicy TwoClassSimulation::heuristicPolicy() const {
	vector<int> bookingLimits(2);
	bookingLimits[0] = (int)floor(min(customers[0].expD * customers[0].consumptionPerReq, (double)C) / customers[0].consumptionPerReq);
	bookingLimits[1] = (int)floor(min(customers[1].expD * customers[1].consumptionPerReq, (double)(C-bookingLimits[0]*customers[0].consumptionPerReq)) / customers[1].consumptionPerReq);
	return bookingLimits;
}

OptionalPolicy TwoClassSimulation::optimalPolicy() const {
	double x = (customers[0].revenuePerReq - customers[1].revenuePerReq) / customers[0].revenuePerReq;
	vector<int> bookingLimits = { C, (int)floor(customers[0].consumptionPerReq * Helpers::invNormal(x, customers[0].expD, customers[0].devD)) };
	return bookingLimits;
}

double MultiClassSimulation::objective(const vector<int>& demands, const vector<int>& bookingLimits) {
	int residualCapacity = C;
	double profit = 0.0;

	for(int ix = numClasses - 1; ix >= 0; ix--) {
		int n = min(demands[ix], (int)floor((bookingLimits[ix] - (C - residualCapacity)) / customers[ix].consumptionPerReq));
		residualCapacity -= n * customers[ix].consumptionPerReq;
		profit += n * customers[ix].revenuePerReq;
	}

	return profit;
}

OptionalPolicy MultiClassSimulation::heuristicPolicy() const {
	vector<int> bookingLimits(numClasses);
	bookingLimits[0] = C;
	int residualCapacity = C - customers[0].expD * customers[0].consumptionPerReq;
	
	for(int j=1; j<numClasses; j++) {
		bookingLimits[j] = (int)floor( min(customers[j].expD * customers[j].consumptionPerReq, (double)residualCapacity) / customers[j].consumptionPerReq );
		residualCapacity -= bookingLimits[j] * customers[j].consumptionPerReq;
	}
	return bookingLimits;
}

OptionalPolicy MultiClassSimulation::optimalPolicy() const {
	return OptionalPolicy();
}

AbstractSimulation::Scenario AbstractSimulation::pickDemands(int scenarioIx, int numScenarios, SamplingType stype) {
    Scenario demands((unsigned long)numClasses);
    int ctr = 0;
    for(Customer &c : customers) {
		double samplingResult = (stype == SamplingType::Descriptive) ? Helpers::pickNormalDescriptive(c.expD, c.devD, scenarioIx, numScenarios) : Helpers::pickNormal(c.expD, c.devD);
		demands[ctr++] = max(0, (int)round(samplingResult));
    }
    return demands;
}

AbstractSimulation::ScenarioList AbstractSimulation::generateScenarios(int ntries, int seed, SamplingType stype) {
	Helpers::resetSeed(seed);
    ScenarioList scenarios((unsigned long) ntries);
    for(int i=0; i<ntries; i++) {
        scenarios[i] = pickDemands(i, ntries, stype);
    }
    return scenarios;
}

Customer::Customer(const json11::Json &obj) :
        name(obj["name"].string_value()),
        expD(obj["expD"].number_value()),
        devD(obj["devD"].number_value()),
        description(obj["description"].string_value()),
        consumptionPerReq(obj["consumptionPerReq"].int_value()),
        revenuePerReq(obj["revenuePerReq"].int_value())
{}
