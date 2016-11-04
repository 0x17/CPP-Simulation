//
// Created by André Schnabel on 30.10.16.
//

#include <algorithm>
#include <cmath>
#include "json11.hpp"

#include "Simulation.h"
#include "Helpers.h"

using namespace std;

std::ostream &operator<<(std::ostream &os, Result const &res) {
	bool first = true;
	os << "profit=" << res.profit << ",";
	os << "bookingLimits=(";
	for (double bl : res.bookingLimits) {
		if (!first)
			os << ", ";
		os << bl;
		first = false;
	}
	os << ")";
	return os;
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
    for(int i=0; i<revenues.size(); i++) {
        revenues[i] = objective(scenarios[i], bookingLimits);
    }
    return revenues;
}

double TwoClassSimulation::objective(const vector<int>& demands, const vector<int>& bookingLimits) {
	int n2 = (int)floor(min(bookingLimits[1], demands[1] * customers[1].consumptionPerReq) / customers[1].consumptionPerReq);
	int n1 = (int)floor(min(demands[0] * customers[0].consumptionPerReq, C - n2 * customers[1].consumptionPerReq));
	return n1 * customers[0].revenuePerReq + n2 * customers[1].revenuePerReq;
}

int TwoClassSimulation::optimalPolicy() {
	double x = (customers[0].revenuePerReq - customers[1].revenuePerReq) / customers[0].revenuePerReq;
	return (int)floor(customers[0].consumptionPerReq * Helpers::invNormal(x, customers[0].expD, customers[0].devD));
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

AbstractSimulation::Scenario AbstractSimulation::pickDemands() {
    Scenario demands((unsigned long)numClasses);
    int ctr = 0;
    for(Customer &c : customers) {
		demands[ctr++] = max(0, (int)round(Helpers::pickNormal(c.expD, c.devD)));
    }
    return demands;
}

AbstractSimulation::ScenarioList AbstractSimulation::generateScenarios(int ntries) {
	Helpers::resetSeed(42);
    ScenarioList scenarios((unsigned long) ntries);
    for(int i=0; i<ntries; i++) {
        scenarios[i] = pickDemands();
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
