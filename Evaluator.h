#pragma once

#include "Simulation.h"

class AbstractEvaluator : public BookingLimitOptimizer
{
public:
	AbstractEvaluator(AbstractSimulation &_sim) : BookingLimitOptimizer("FullEnumeration", _sim) {}
	virtual ~AbstractEvaluator() {}

	virtual ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) const = 0;
	virtual Result computeOptimum(AbstractSimulation::ScenarioList &scenarios) const = 0;

	static Result extractOptimumFromList(const ResultList &results, bool printOpts = false);

	Result solve(AbstractSimulation::ScenarioList& scenarios) override;
};

class Evaluator2D : public AbstractEvaluator {
public:
	explicit Evaluator2D(AbstractSimulation &_sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) const override;
	Result computeOptimum(AbstractSimulation::ScenarioList& scenarios) const override;
};

class Evaluator3D : public AbstractEvaluator {
public:
	explicit Evaluator3D(AbstractSimulation &_sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) const override;
	Result computeOptimum(AbstractSimulation::ScenarioList& scenarios) const override;
};

class EvaluatorMultiDimensional : public AbstractEvaluator {
public:
	explicit EvaluatorMultiDimensional(AbstractSimulation &_sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(AbstractSimulation::ScenarioList &scenarios) const override;
	Result computeOptimum(AbstractSimulation::ScenarioList& scenarios) const override;
};

