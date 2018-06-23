#pragma once

#include <list>
#include <gurobi_c++.h>

#include "Simulation.h"

struct DisplayVar {
	std::string name;
	GRBVar v;
};

class GurobiOptimizer : public BookingLimitOptimizer {
public:
	explicit GurobiOptimizer(const AbstractSimulation& _sim) : BookingLimitOptimizer("Gurobi", _sim) {}
	Result solve(const DemandScenarioList& scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) override;

	void addDisplayVar(const DisplayVar& v);
	void addDisplayVars(const std::vector<DisplayVar>& vars);
	void addDisplayVars(const std::vector<std::string> &names, const std::vector<GRBVar>& grbVars);
	void addDisplayVarsContiguous(const std::string &prefix, const std::vector<GRBVar>& grbVars);

private:
	void printDisplayVarResults();
	template<class Func>
	Result solveCommon(GurobiOptimizer &solver, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc, Func modelBuilder);

	std::list<DisplayVar> displayVars;
};