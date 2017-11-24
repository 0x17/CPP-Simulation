//
// Created by André Schnabel on 30.10.16.
//

#ifndef CPP_SIMULATION_SIMULATION_H
#define CPP_SIMULATION_SIMULATION_H

#include <boost/optional.hpp>
#include "json11.hpp"
#include "Matrix.h"

struct Customer {
    std::string name;
    double expD, devD;
    std::string description;
    double consumptionPerReq, revenuePerReq;

    explicit Customer(const json11::Json &obj);
};

using OptionalPolicy = boost::optional<std::vector<int>>;

using Scenario = std::vector<int>;
using ScenarioList = Matrix<int>;

class AbstractSimulation {
public:
	using LUTList = std::vector<std::vector<double>>;

    explicit AbstractSimulation(const std::string &dataFilename);
	virtual ~AbstractSimulation() {}

	enum class SamplingType {
		Random,
		Descriptive
	};

    Scenario pickDemands(int scenarioIx, int numScenarios);
	Scenario pickDemandsDescriptive(int scenarioIx, int numScenarios, LUTList& lutList);
    ScenarioList generateScenarios(int ntries, int seed, SamplingType stype = SamplingType::Descriptive);
    std::vector<double> runSimulation(const std::vector<int>& bookingLimits, const ScenarioList &scenarios) const;
	double averageRevenueOfSimulation(const std::vector<int>& bookingLimits, const ScenarioList& scenarios) const;

	virtual double objective(const std::vector<int> &demands, const std::vector<int> &bookingLimits) const = 0;

	virtual OptionalPolicy optimalPolicy() const = 0;
	virtual OptionalPolicy heuristicPolicy() const = 0;

	int getC() const { return C; }
	int getNumClasses() const { return (int)customers.size();  }

	Customer getCustomer(int ix) const { return customers[ix];  }

	static std::vector<double> statisticalMeansOfScenarios(ScenarioList &scenarios);
	static std::vector<double> statisticalStandardDeviationsOfScenarios(ScenarioList &scenarios);

protected:
	int C;
    std::vector<Customer> customers;
    int numClasses;
};

struct Result {
	std::vector<int> bookingLimits;
	double profit;

	Result() : bookingLimits(0), profit(0) {}
    explicit Result(int numClasses) : bookingLimits((unsigned long)numClasses), profit(0) {}
	Result(const std::vector<int>& booking_limits, double profit) : bookingLimits(booking_limits), profit(profit) {}

	std::string toString() const;
};

using ResultList = std::vector<Result>;

class BookingLimitOptimizer {
	const bool useHeuristicStart = false;
public:
	BookingLimitOptimizer(std::string _name, const AbstractSimulation &_sim) : sim(_sim), name(_name) {
		if (useHeuristicStart) {
			heuristicBookingLimits = sim.heuristicPolicy();
		}
	}
	virtual ~BookingLimitOptimizer() {}
	virtual Result solve(const ScenarioList& scenarios) = 0;

	std::string getName() const { return name; }

protected:
	const AbstractSimulation &sim;
	OptionalPolicy heuristicBookingLimits;
private:
	std::string name;
};

class TwoClassSimulation : public AbstractSimulation {
public:
    explicit TwoClassSimulation(const std::string &dataFilename) : AbstractSimulation(dataFilename) {}
    double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) const override;
	OptionalPolicy heuristicPolicy() const override;
	OptionalPolicy optimalPolicy() const override;
};

class MultiClassSimulation : public AbstractSimulation {
public:
    explicit MultiClassSimulation(const std::string &dataFilename = "multi_data.json") : AbstractSimulation(dataFilename) {}
    double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) const override;
	OptionalPolicy heuristicPolicy() const override;
	OptionalPolicy optimalPolicy() const override;
	double eosConsumption(int j, int u) const;
};


#endif //CPP_SIMULATION_SIMULATION_H
