#pragma once

#include "Simulation.h"

class AbstractEvaluator : public BookingLimitOptimizer
{
public:
	explicit AbstractEvaluator(const AbstractSimulation &_sim) : BookingLimitOptimizer("FullEnumeration", _sim) {}
	virtual ~AbstractEvaluator() = default;

	virtual ResultList collectResults(const DemandScenarioList &scenarios) const = 0;
	virtual Result computeOptimum(const DemandScenarioList &scenarios) const = 0;

	static Result extractOptimumFromList(const ResultList &results, bool printOpts = false);

	Result solve(const DemandScenarioList& scenarios) override;
};

class Evaluator2D : public AbstractEvaluator {
public:
	explicit Evaluator2D(const AbstractSimulation &_sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(const DemandScenarioList &scenarios) const override;
	Result computeOptimum(const DemandScenarioList& scenarios) const override;
};

class Evaluator3D : public AbstractEvaluator {
public:
	explicit Evaluator3D(const AbstractSimulation &_sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(const DemandScenarioList &scenarios) const override;
	Result computeOptimum(const DemandScenarioList& scenarios) const override;
};

class EvaluatorMultiDimensional : public AbstractEvaluator {
public:
	explicit EvaluatorMultiDimensional(const AbstractSimulation &_sim) : AbstractEvaluator(_sim) {}
	ResultList collectResults(const DemandScenarioList &scenarios) const override;
	Result computeOptimum(const DemandScenarioList& scenarios) const override;
};

