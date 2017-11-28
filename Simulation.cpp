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
#include <boost/algorithm/clamp.hpp>

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
    auto obj = json11::Json::parse(Helpers::slurp(dataFilename), errMsg);
	if(!errMsg.empty()) {
		throw runtime_error("Unable to parse " + dataFilename + " error = " + errMsg);
	}
    C = obj["capacity"].int_value();
	for(const auto &c : obj["clients"].array_items()) {
		customers.emplace_back(c);
	}
    numClasses = (int)customers.size();
}

vector<double> AbstractSimulation::runSimulation(const vector<int> &bookingLimits, const DemandScenarioList &scenarios) const {
    vector<double> revenues(static_cast<unsigned long>(scenarios.getM()));
    for(int i=0; i<revenues.size(); i++) {
        revenues[i] = objective(scenarios.row(i), bookingLimits);
    }
    return revenues;
}

double AbstractSimulation::averageRevenueOfSimulation(const std::vector<int>& bookingLimits, const DemandScenarioList& scenarios) const {
	return Helpers::vecAverage(runSimulation(bookingLimits, scenarios));
}

double AbstractSimulation::conditionalValueAtRiskOfSimulationResult(double alpha, const std::vector<double> &revenues) const {
	int numWorstScenarios = round((1 - alpha) * revenues.size());
	std::vector<double> rcopy(revenues);
	sort(rcopy.begin(), rcopy.end());
	return Helpers::vecAverageSubRange(rcopy, 0, numWorstScenarios);
}

double AbstractSimulation::weightedProfitAndCVaRofSimulationResult(double profitWeight, double alpha, const std::vector<double>& revenues) const {
	double cvar = conditionalValueAtRiskOfSimulationResult(alpha, revenues);
	double avgProfit = Helpers::vecAverage(revenues);
	return profitWeight * avgProfit + (1 - profitWeight) * cvar;
}

double AbstractSimulation::objectiveWithGlobalSettings(const std::vector<int>& bookingLimits, const DemandScenarioList& scenarios) const {
	auto revenues = runSimulation(bookingLimits, scenarios);
	return globals::CONDITIONAL_VALUE_AT_RISK_ENABLED ? weightedProfitAndCVaRofSimulationResult(globals::PROFIT_WEIGHT, globals::ALPHA, revenues) : Helpers::vecAverage(revenues);
}

vector<double> AbstractSimulation::statisticalMeansOfScenarios(DemandScenarioList& scenarios) {
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

vector<double> AbstractSimulation::statisticalStandardDeviationsOfScenarios(DemandScenarioList& scenarios) {
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
	int n2 = (int)floor(min((double)bookingLimits[1], demands[1] * customers[1].consumptionPerReqMean) / customers[1].consumptionPerReqMean);
	int n1 = (int)floor(min(demands[0] * customers[0].consumptionPerReqMean, C - n2 * customers[1].consumptionPerReqMean));
	return n1 * customers[0].revenuePerReq + n2 * customers[1].revenuePerReq;
}

OptionalPolicy TwoClassSimulation::heuristicPolicy() const {
	vector<int> bookingLimits(2);
	bookingLimits[0] = (int)floor(min(customers[0].expD * customers[0].consumptionPerReqMean, (double)C) / customers[0].consumptionPerReqMean);
	bookingLimits[1] = (int)floor(min(customers[1].expD * customers[1].consumptionPerReqMean, (C-bookingLimits[0]*customers[0].consumptionPerReqMean)) / customers[1].consumptionPerReqMean);
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
			n = min(demands[ix], (int)floor((bookingLimits[ix] - (C - residualCapacity)) / customers[ix].consumptionPerReqMean));
			residualCapacity -= n * customers[ix].consumptionPerReqMean;
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
		bookingLimits[j] = (int)floor( min(customers[j].expD * customers[j].consumptionPerReqMean, (double)bookingLimits[j-1]) );
	}
	return bookingLimits;
}

OptionalPolicy MultiClassSimulation::optimalPolicy() const {
	return OptionalPolicy();
}

double MultiClassSimulation::eosConsumption(int j, int u) const {
	double deltaX = customers[j].expD;
	double deltaY = customers[j].consumptionPerReqMean * 0.2;
	return max(1.0, customers[j].consumptionPerReqMean - (boost::math::erf(4 * u / deltaX - 2) * deltaY / 2.0 + deltaY / 2.0));
}

template<class T>
using Scenario = std::vector<T>;

template<class T>
using ScenarioList = Matrix<T>;

template<class T>
Scenario<T> pickScenario(int nclasses, const vector<DistParameters> &distParams, bool doRound = false) {
	Scenario<T> scenario(nclasses);
	for (int classIx = 0; classIx < nclasses; classIx++) {
		double samplingResult = Helpers::pickNormal(distParams[classIx].mean, distParams[classIx].stddev);
		scenario[classIx] = doRound ? max(0, (int)round(samplingResult)) : samplingResult;
	}
	return scenario;
}

template<class T>
Scenario<T> pickScenarioDescriptive(int nclasses, LUTList &lutList, bool doRound = false) {
	Scenario<T> scenario(nclasses);
	for (int classIx = 0; classIx < nclasses; classIx++) {
		double samplingResult = Helpers::pickNextWithLUT(lutList[classIx]);
		scenario[classIx] = doRound ? max(0, (int)round(samplingResult)) : samplingResult;
	}
	return scenario;
}

template<class T>
ScenarioList<T> generateScenarios(int nclasses, const vector<DistParameters> &distParams, int ntries, int seed, SamplingType stype, bool doRound = false) {
	Helpers::resetSeed(seed);

	ScenarioList<T> scenarios(ntries, static_cast<int>(nclasses));

	if (stype == SamplingType::Descriptive) {
		LUTList lookupTables = generateLookupTableList(nclasses, ntries, distParams);

		for (int i = 0; i<ntries; i++) {
			scenarios.setRow(i, pickScenarioDescriptive<T>(nclasses, lookupTables, doRound));
		}
	}
	else {
		for (int i = 0; i<ntries; i++) {
			scenarios.setRow(i, pickScenario<T>(nclasses, distParams, doRound));
		}
	}

	return scenarios;
}

DemandScenario AbstractSimulation::pickDemands() const {
	return pickScenario<int>(numClasses, demandDistributionParametersForCustomers(), true);
}

DemandScenario AbstractSimulation::pickDemandsDescriptive(LUTList &lutList) const {
	return pickScenarioDescriptive<int>(customers.size(), lutList, true);
}

DemandScenarioList AbstractSimulation::generateDemandScenarios(int ntries, int seed, SamplingType stype) const {
	return generateScenarios<int>(customers.size(), demandDistributionParametersForCustomers(), ntries, seed, stype, true);
}

ConsumptionScenario AbstractSimulation::pickConsumptions() const {
	return pickScenario<double>(numClasses, consumptionDistributionParametersForCustomers(), false);
}

ConsumptionScenario AbstractSimulation::pickConsumptionsDescriptive(LUTList& lutList) const {
	return pickScenarioDescriptive<double>(customers.size(), lutList, false);
}

ConsumptionScenarioList AbstractSimulation::generateConsumptionScenarios(int ntries, int seed, SamplingType stype) const {
	return generateScenarios<double>(customers.size(), consumptionDistributionParametersForCustomers(), ntries, seed, stype, false);
}

std::vector<DistParameters> AbstractSimulation::demandDistributionParametersForCustomers() const {
    return Helpers::constructVector<DistParameters>(getNumClasses(), [this](int ix) {
        return DistParameters { customers[ix].expD, customers[ix].devD };
    });
}

std::vector<DistParameters> AbstractSimulation::consumptionDistributionParametersForCustomers() const {
    return Helpers::constructVector<DistParameters>(getNumClasses(), [this](int ix) {
        return DistParameters { customers[ix].consumptionPerReqMean, customers[ix].consumptionPerReqStdDev };
    });
}

LUTList generateLookupTableList(int nclasses, int ntries, const std::vector<DistParameters> &distParams) {
    LUTList lookupTables(nclasses);
    for (int i = 0; i < nclasses; i++) {
        lookupTables[i] = Helpers::generateNormalDistributionDescriptiveSamplingLUT(ntries, distParams[i].mean, distParams[i].stddev);
    }
    return lookupTables;
}

Customer::Customer(const json11::Json &obj) :
        name(obj["name"].string_value()),
        expD(obj["expD"].number_value()),
        devD(obj["devD"].number_value()),
        description(obj["description"].string_value()),
        consumptionPerReqMean(obj["consumptionPerReqMean"].number_value()),
		consumptionPerReqStdDev(obj["consumptionPerReqStdDev"].is_number() ? obj["consumptionPerReqStdDev"].number_value() : 0.0),
        revenuePerReq(obj["revenuePerReq"].number_value())
{}
