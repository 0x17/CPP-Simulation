//
// Created by André Schnabel on 30.10.16.
//

#pragma once

#include <boost/optional.hpp>
#include "json11.hpp"
#include "Matrix.h"
#include "Helpers.h"
#include "Matrix3D.h"

struct Toggles {
	bool economyOfScale, CVaR, stochasticConsumptions;

	explicit Toggles(const std::string& filename = "Toggles.json");

	std::vector<bool> toVec() const;
};

struct Customer {
    std::string name;
    double expD, devD;
    std::string description;
	double consumptionPerReqMean, consumptionPerReqStdDev;
	double revenuePerReq;

    explicit Customer(const json11::Json &obj);
};

using OptionalPolicy = boost::optional<std::vector<int>>;

using DemandScenario = std::vector<int>;
using DemandScenarioList = Matrix<int>;

using ConsumptionScenario = std::vector<double>;
using ConsumptionScenarioList = Matrix<double>;
using ConsumptionScenarioFunc = Matrix3D<double>;

using LUTList = std::vector<std::vector<double>>;

enum class SamplingType {
	Random,
	Descriptive
};

struct DistParameters {
    double mean, stddev;
};

LUTList generateLookupTableList(int nclasses, int ntries, const std::vector<DistParameters> &distParams);

class AbstractSimulation {
public:
    explicit AbstractSimulation(const std::string &dataFilename, Toggles _toggles);
	virtual ~AbstractSimulation() = default;

    DemandScenario pickDemands() const;
	DemandScenario pickDemandsDescriptive(LUTList& lutList) const;
    DemandScenarioList generateDemandScenarios(int ntries, int seed, SamplingType stype = SamplingType::Descriptive) const;

	ConsumptionScenario pickConsumptions(int u) const;
	ConsumptionScenario pickConsumptionsDescriptive(LUTList& lutList) const;
	ConsumptionScenarioList generateConsumptionScenarios(int u, int ntries, int seed, SamplingType stype = SamplingType::Descriptive) const;
	ConsumptionScenarioFunc generateConsumptionScenarioFunc(int ntries, int seed, SamplingType stype = SamplingType::Descriptive, int maxAccept = -1) const;

    std::vector<double> runSimulation(const std::vector<int>& bookingLimits, const DemandScenarioList &scenarios) const;

	double averageRevenueOfSimulation(const std::vector<int>& bookingLimits, const DemandScenarioList& scenarios) const;
	double conditionalValueAtRiskOfSimulationResult(double alpha, const std::vector<double> &revenues) const;
	double weightedProfitAndCVaRofSimulationResult(double profitWeight, double alpha, const std::vector<double> &revenues) const;

	double objectiveWithCVarOption(const std::vector<int>& bookingLimits, const DemandScenarioList& scenarios) const;

	virtual double objective(const std::vector<int> &demands, const std::vector<int> &bookingLimits) const = 0;

	virtual OptionalPolicy optimalPolicy() const = 0;
	virtual OptionalPolicy heuristicPolicy() const = 0;

	int getC() const { return C; }
	int getNumClasses() const { return (int)customers.size();  }

	Customer getCustomer(int ix) const { return customers[ix];  }

    std::vector<DistParameters> demandDistributionParametersForCustomers() const;
    std::vector<DistParameters> consumptionDistributionParametersForCustomers(int u) const;

	static std::vector<double> statisticalMeansOfScenarios(DemandScenarioList &scenarios);
	static std::vector<double> statisticalStandardDeviationsOfScenarios(DemandScenarioList &scenarios);

	virtual DistParameters eosConsumptionDistributionParametersForCustomer(int j, int u) const = 0;

	Toggles getToggles() const;

	void setRiskAversionParameters(double _alpha, double _psi);
	double getAlpha() const { return alpha;  }
	double getPsi() const { return psi; }

protected:
	Toggles toggles;
	int C;
    std::vector<Customer> customers;
    int numClasses;
	double alpha, psi;
};

struct Result {
	std::vector<int> bookingLimits;
	double profit;

	Result();
	explicit Result(int numClasses);
	Result(const std::vector<int>& booking_limits, double profit);

	std::string toString() const;
};

using ResultList = std::vector<Result>;

class BookingLimitOptimizer {
	const bool useHeuristicStart = false;
public:
	BookingLimitOptimizer(std::string _name, const AbstractSimulation& _sim);
	virtual ~BookingLimitOptimizer() = default;

	virtual Result solve(const DemandScenarioList& scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) = 0;

	std::string getName() const;

protected:
	const AbstractSimulation &sim;
	OptionalPolicy heuristicBookingLimits;
private:
	std::string name;
};

class TwoClassSimulation : public AbstractSimulation {
public:
	explicit TwoClassSimulation(const std::string& dataFilename, Toggles _toggles);
    double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) const override;
	OptionalPolicy heuristicPolicy() const override;
	OptionalPolicy optimalPolicy() const override;
	DistParameters eosConsumptionDistributionParametersForCustomer(int j, int u) const override;
};

class MultiClassSimulation : public AbstractSimulation {
public:
	explicit MultiClassSimulation(const std::string& dataFilename, Toggles _toggles);
    double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) const override;
	OptionalPolicy heuristicPolicy() const override;
	OptionalPolicy optimalPolicy() const override;
	double eosConsumption(int j, int u) const;
	DistParameters eosConsumptionDistributionParametersForCustomer(int j, int u) const override;
};

