//
// Created by André Schnabel on 30.10.16.
//

#ifndef CPP_SIMULATION_SIMULATION_H
#define CPP_SIMULATION_SIMULATION_H

#include <list>
#include <boost/optional.hpp>
#include "json11.hpp"

struct Customer {
    std::string name;
    double expD, devD;
    std::string description;
    int consumptionPerReq, revenuePerReq;

    Customer(const json11::Json &obj);
};

using OptionalPolicy = boost::optional<std::vector<int>>;

class AbstractSimulation {
public:
	using Scenario = std::vector<int>;
	using ScenarioList = std::vector<Scenario>;

	AbstractSimulation(const std::string &dataFilename);
	virtual ~AbstractSimulation() {}

	enum class SamplingType {
		Random,
		Descriptive
	};

    Scenario pickDemands(int scenarioIx, int numScenarios, SamplingType stype = SamplingType::Descriptive);
    ScenarioList generateScenarios(int ntries, int seed, SamplingType stype = SamplingType::Descriptive);
    std::vector<double> runSimulation(const std::vector<int> &bookingLimits, ScenarioList &scenarios);
	double averageRevenueOfSimulation(const std::vector<int>& bookingLimits, ScenarioList& scenarios);

	virtual double objective(const std::vector<int> &demands, const std::vector<int> &bookingLimits) = 0;

	virtual OptionalPolicy optimalPolicy() const = 0;
	virtual OptionalPolicy heuristicPolicy() const = 0;

	int getC() const { return C; }
	int getNumClasses() const { return (int)customers.size();  }

	Customer &getCustomer(int ix) { return customers[ix];  }

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
	Result(int nclasses) : bookingLimits(nclasses), profit(0) {}
	Result(const std::vector<int>& booking_limits, double profit) : bookingLimits(booking_limits), profit(profit) {}

	std::string toString() const;
};

using ResultList = std::vector<Result>;

class BookingLimitOptimizer {
	const bool useHeuristicStart = true;
public:
	BookingLimitOptimizer(std::string _name, AbstractSimulation &_sim) : sim(_sim), name(_name) {
		if (useHeuristicStart) {
			heuristicBookingLimits = sim.heuristicPolicy();
		}
	}
	virtual ~BookingLimitOptimizer() {}
	virtual Result solve(std::vector<std::vector<int>>& scenarios) = 0;

	std::string getName() const { return name; }

protected:
	AbstractSimulation &sim;
	OptionalPolicy heuristicBookingLimits;
private:
	std::string name;
};

class TwoClassSimulation : public AbstractSimulation {
public:
    TwoClassSimulation(const std::string &dataFilename) : AbstractSimulation(dataFilename) {}
	virtual double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) override;
	OptionalPolicy heuristicPolicy() const override;
	OptionalPolicy optimalPolicy() const override;
};

class MultiClassSimulation : public AbstractSimulation {
public:
    MultiClassSimulation(const std::string &dataFilename = "multi_data.json") : AbstractSimulation(dataFilename) {}
	virtual double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) override;
	OptionalPolicy heuristicPolicy() const override;
	OptionalPolicy optimalPolicy() const override;
};


#endif //CPP_SIMULATION_SIMULATION_H
