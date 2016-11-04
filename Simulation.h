//
// Created by André Schnabel on 30.10.16.
//

#ifndef CPP_SIMULATION_SIMULATION_H
#define CPP_SIMULATION_SIMULATION_H

#include <list>
#include "json11.hpp"

struct Customer {
    std::string name;
    double expD, devD;
    std::string description;
    int consumptionPerReq, revenuePerReq;

    Customer(const json11::Json &obj);
};

class AbstractSimulation {
public:
	using Scenario = std::vector<int>;
	using ScenarioList = std::vector<Scenario>;

	AbstractSimulation(const std::string &dataFilename);
	virtual ~AbstractSimulation() {}

    Scenario pickDemands();
    ScenarioList generateScenarios(int ntries);
    std::vector<double> runSimulation(const std::vector<int> &bookingLimits, ScenarioList &scenarios);

    virtual double objective(const std::vector<int> &demands, const std::vector<int> &bookingLimits) = 0;

	int getC() const { return C; }
	int getNumClasses() const { return (int)customers.size();  }

	Customer &getCustomer(int ix) { return customers[ix];  }

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
};

using ResultList = std::vector<Result>;

std::ostream &operator<<(std::ostream &os, Result const &res);

class BookingLimitOptimizer {
public:
	BookingLimitOptimizer(std::string _name, AbstractSimulation &_sim) : sim(_sim), name(_name) {}
	virtual ~BookingLimitOptimizer() {}
	virtual Result solve(std::vector<std::vector<int>>& scenarios) = 0;

	std::string getName() const { return name; }

protected:
	AbstractSimulation &sim;
private:
	std::string name;
};

class TwoClassSimulation : public AbstractSimulation {
public:
    TwoClassSimulation(const std::string &dataFilename) : AbstractSimulation(dataFilename) {}
	virtual double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) override;
	int optimalPolicy();
};

class MultiClassSimulation : public AbstractSimulation {
public:
    MultiClassSimulation() : AbstractSimulation("multi_data.json") {}
	virtual double objective(const std::vector<int>& demands, const std::vector<int>& bookingLimits) override;
};


#endif //CPP_SIMULATION_SIMULATION_H
