#pragma once
#include <tuple>
#include <vector>
#include "Simulation.h"

class AbstractSimulation;

class AbstractEvaluator
{
public:
	struct Result {
		std::vector<int> bookingLimits;
		double profit;

		Result() : bookingLimits(0), profit(0) {}
		Result(int nclasses) : bookingLimits(nclasses), profit(0) {}
	};
	using ResultList = std::vector<Result>;

	AbstractEvaluator(AbstractSimulation &_sim) : sim(_sim) {}
	virtual ~AbstractEvaluator() {}

	virtual ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) = 0;

	static Result computeOpt(const ResultList &results, bool printOpts = false);

protected:
	AbstractSimulation &sim;
};

std::ostream &operator<<(std::ostream &os, AbstractEvaluator::Result const &res);

class Evaluator2D : public AbstractEvaluator {
public:
	explicit Evaluator2D(TwoClassSimulation& _sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) override;
};

class Evaluator3D : public AbstractEvaluator {

public:
	explicit Evaluator3D(MultiClassSimulation& _sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) override;
};

