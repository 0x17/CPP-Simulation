#pragma once
#include <tuple>
#include <vector>
#include "Simulation.h"

class AbstractSimulation;

class AbstractEvaluator : public BookingLimitOptimizer
{
public:
	AbstractEvaluator(AbstractSimulation &_sim) : BookingLimitOptimizer("FullEnumeration", _sim) {}
	virtual ~AbstractEvaluator() {}

	virtual ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) = 0;

	static Result computeOpt(const ResultList &results, bool printOpts = false);

	virtual Result solve(AbstractSimulation::ScenarioList& scenarios) override;
};

class Evaluator2D : public AbstractEvaluator {
public:
	explicit Evaluator2D(AbstractSimulation &_sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) override;
};

class EvaluatorMultiDimensional : public AbstractEvaluator {

public:
	explicit EvaluatorMultiDimensional(AbstractSimulation &_sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) override;
};

