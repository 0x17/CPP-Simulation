#include <gurobi_c++.h>

#include "GurobiSolver.h"
#include "Simulation.h"
#include "Matrix.h"
#include "Helpers.h"

using namespace std;

Result GurobiOptimizer::solve(std::vector<std::vector<int>>& scenarios) {
	Result res;

	int J = sim.getNumClasses(),
		C = sim.getC(),
		S = (int)scenarios.size();

	auto cj = [&](int j) { return sim.getCustomer(j).consumptionPerReq; };
	auto rj = [&](int j) { return sim.getCustomer(j).revenuePerReq; };

	GRBEnv env;
	env.set(GRB_DoubleParam_MIPGap, 0.0);
	env.set(GRB_DoubleParam_TimeLimit, 5);

	GRBModel model(env);

	vector<GRBVar> bcj = Helpers::constructVector<GRBVar>(J, [&](int j) { return model.addVar(j == 0 ? C : 0, C, 0.0, GRB_INTEGER, "bcj" + to_string(j)); });
	Matrix<GRBVar> njs(J, S, [&](int j, int s) { return model.addVar(0, floor((double)C / (double)cj(j)), 0.0, GRB_CONTINUOUS, "njs" + to_string(j) + "," + to_string(s)); });
	Matrix<GRBVar> rcapjs(J, S, [&](int j, int s) { return model.addVar(0, C, 0.0, GRB_CONTINUOUS, "rcapjs" + to_string(j) + "," + to_string(s)); });
	Matrix<GRBVar> kjs(J, S, [&](int j, int s) { return model.addVar(0, floor((double)C / (double)cj(j)), 0.0, GRB_INTEGER, "kjs" + to_string(j) + "," + to_string(s)); });

	GRBLinExpr revSum = 0.0;
	for (int j = 0; j < J; j++)
		for (int s = 0; s < S; s++)
			revSum += njs(j, s) * (double)rj(j);

	model.setObjective(1.0 / (double)S * revSum, GRB_MAXIMIZE);

	for (int j = 0; j < J; j++) {
		if (j + 1 < J)
			model.addConstr(bcj[j] >= bcj[j + 1]);

		for (int s = 0; s < S; s++) {
			model.addConstr(kjs(j, s) * cj(j) <= rcapjs(j, s));
			model.addConstr(kjs(j, s) * cj(j) >= rcapjs(j, s) - cj(j) + 1);

			GRBLinExpr previouslyUtilizedCap = 0.0;
			for (int i = j + 1; i < J; i++)
				previouslyUtilizedCap += njs(i, s) * cj(i);
			model.addConstr(rcapjs(j, s) == (bcj[j] - previouslyUtilizedCap));

			GRBVar kjsarr[] = { kjs(j,s) };
			model.addGenConstrMin(njs(j, s), kjsarr, 1, scenarios[s][j]);
		}
	}

	model.update();

	try {
		model.optimize();

		res.profit = model.get(GRB_DoubleAttr_ObjVal);
		res.bookingLimits = Helpers::constructVector<int>(J, [&](int j) { return (int)bcj[j].get(GRB_DoubleAttr_X); });;

		cout << "Optimality: " << (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) << endl;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}

	return res;
}
