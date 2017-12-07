#include <gurobi_c++.h>
#include <cmath>
#include <boost/math/special_functions.hpp>

#include "GurobiSolver.h"
#include "Simulation.h"
#include "Matrix.h"
#include "Helpers.h"
#include "Globals.h"
#include <numeric>

using namespace std;

class CustomCallback : public GRBCallback {
public:
	CustomCallback();
private:
	void callback() override;
	Helpers::Tracer tr;
	Stopwatch sw;
};

CustomCallback::CustomCallback() : tr("GurobiTrace") {
	sw.start();
}

void CustomCallback::callback() {
	if (where == GRB_CB_MIP) {
		tr.trace(/*getDoubleInfo(GRB_CB_RUNTIME)*/ sw.look(), static_cast<float>(getDoubleInfo(GRB_CB_MIP_OBJBST)));
	}
}

struct BookingAcceptVars {
	vector<GRBVar> bcj;
	Matrix<GRBVar> njs;
};

BookingAcceptVars modelBuilderForOldFormulation(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc);
BookingAcceptVars modelBuilderForNewFormulation(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc);
BookingAcceptVars modelBuilderForEconomiesOfScale(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc);
BookingAcceptVars modelBuilderForStochasticConsumptions(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc);
BookingAcceptVars modelBuilderForStochasticEconomiesOfScale(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc);
BookingAcceptVars modelBuilderForConditionalValueAtRisk(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc);

Result GurobiOptimizer::solve(const DemandScenarioList& scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) {
	auto assertSingleToggle = [](const Toggles &toggleObj) {
		const vector<bool> &toggles = toggleObj.toVec();
		auto b2int = [](bool b) -> int { return b ? 1 : 0;  };
		vector<int> nums(toggles.size());
		std::transform(toggles.begin(), toggles.end(), nums.begin(), b2int);
		int ntoggles = std::accumulate(nums.begin(), nums.end(), 0);
		if (ntoggles > 1) {
			throw runtime_error("More than one toggle!");
		}
	};

	const Toggles toggles = sim.getToggles();

	assertSingleToggle(toggles);

	auto chooseModelBuilder = [](const Toggles &toggles) {
		if(toggles.economyOfScale && toggles.stochasticConsumptions) return modelBuilderForStochasticEconomiesOfScale;
		if(toggles.economyOfScale) return modelBuilderForEconomiesOfScale;
		//if(toggles.conditionalValueAtRisk) return modelBuilderForConditionalValueAtRisk;
		if(toggles.stochasticConsumptions) return modelBuilderForStochasticConsumptions;
		return modelBuilderForNewFormulation;
	};

	auto modelBuilder = chooseModelBuilder(toggles);
	return solveCommon(scenarios, consumptionScenarioFunc, modelBuilder);
}

void addConditionalValueAtRiskConsiderationToModel(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const Matrix<GRBVar> &njs);

template<class Func>
Result GurobiOptimizer::solveCommon(const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc, Func modelBuilder) {
	GRBEnv env;
	env.set(GRB_DoubleParam_MIPGap, 0.0);
	env.set(GRB_DoubleParam_TimeLimit, /*GRB_INFINITY*/ globals::TIME_LIMIT);
	env.set(GRB_IntParam_Threads, 1);

	GRBModel model(env);

	BookingAcceptVars pair = modelBuilder(sim, model, scenarios, consumptionScenarioFunc);
	auto bcj = pair.bcj;
    auto njs = pair.njs;

    if(sim.getToggles().CVaR) {
        addConditionalValueAtRiskConsiderationToModel(sim, model, scenarios, njs);
    }

	if (heuristicBookingLimits) {
		for (int j = 0; j<scenarios.getN(); j++) {
			bcj[j].set(GRB_DoubleAttr_Start, (*heuristicBookingLimits)[j]);
		}
	}

	CustomCallback callback;

	model.update();
	model.setCallback(&callback);

	Result res;

	try {
		model.optimize();

		res.profit = model.get(GRB_DoubleAttr_ObjVal);
		res.bookingLimits = Helpers::constructVector<int>(scenarios.getN(), [&](int j) { return (int)round(bcj[j].get(GRB_DoubleAttr_X)); });;

		cout << "Optimality: " << (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) << endl;

		double secondsElapsed = model.get(GRB_DoubleAttr_Runtime);
		//Helpers::spitAppend(to_string(S) + ";" + to_string(secondsElapsed) + "\n", "solvetimeforntries.txt");
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}

	return res;
}

template<class Func>
GRBLinExpr sum(int ub, Func f) {
	GRBLinExpr result = 0.0;
	for(int i=0; i < ub; i++) {
		result += f(i);
	}
	return result;
}

template<class Func>
GRBLinExpr sum(int lb, int ub, Func f) {
	GRBLinExpr result = 0.0;
	for(int i=lb; i < ub; i++) {
		result += f(i);
	}
	return result;
}

template<class Func>
GRBLinExpr sum2D(int ub1, int ub2, Func f) {
	GRBLinExpr result = 0.0;
	for(int i=0; i < ub1; i++) {
		for(int j=0; j < ub2; j++) {
			result += f(i,j);
		}
	}
	return result;
}

template<class Func>
GRBLinExpr sum2D(int lb1, int ub1, int lb2, int ub2, Func f) {
	GRBLinExpr result = 0.0;
	for(int i=lb1; i < ub1; i++) {
		for(int j=lb2; j < ub2; j++) {
			result += f(i,j);
		}
	}
	return result;
}

template<class Func>
GRBLinExpr sum3D(int ub1, int ub2, int ub3, Func f) {
	GRBLinExpr result = 0.0;
	for(int i=0; i < ub1; i++) {
		for(int j=0; j < ub2; j++) {
			for(int k=0; k < ub3; k++) {
				result += f(i, j, k);
			}
		}
	}
	return result;
}

BookingAcceptVars modelBuilderForOldFormulation(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) {
	int J = sim.getNumClasses(),
		C = sim.getC(),
		S = scenarios.getM();

	auto cj = [&](int j) { return sim.getCustomer(j).consumptionPerReqMean; };
	auto rj = [&](int j) { return sim.getCustomer(j).revenuePerReq; };

	vector<GRBVar> bcj = Helpers::constructVector<GRBVar>(J, [&](int j) { return model.addVar(j == 0 ? C : 0, C, 0.0, GRB_INTEGER, "bcj" + to_string(j)); });
	Matrix<GRBVar> njs(J, S, [&](int j, int s) { return model.addVar(0, std::floor((double)C / (double)cj(j)), 0.0, GRB_CONTINUOUS, "njs" + to_string(j) + "," + to_string(s)); });
	Matrix<GRBVar> rcapjs(J, S, [&](int j, int s) { return model.addVar(0, C, 0.0, GRB_CONTINUOUS, "rcapjs" + to_string(j) + "," + to_string(s)); });
	Matrix<GRBVar> kjs(J, S, [&](int j, int s) { return model.addVar(0, std::floor((double)C / (double)cj(j)), 0.0, GRB_INTEGER, "kjs" + to_string(j) + "," + to_string(s)); });

	model.setObjective(1.0 / (double)S * sum2D(J, S, [&](int j, int s) {  return njs(j, s) * (double)rj(j); }), GRB_MAXIMIZE);

	for (int j = 0; j < J; j++) {
		if (j + 1 < J)
			model.addConstr(bcj[j] >= bcj[j + 1]);

		for (int s = 0; s < S; s++) {
			model.addConstr(kjs(j, s) * cj(j) <= rcapjs(j, s));
			model.addConstr((kjs(j, s) + 1) * cj(j) >= rcapjs(j, s) + globals::EPSILON);
			//model.addConstr(kjs(j, s) * cj(j) >= rcapjs(j, s) - cj(j) + globals::EPSILON);

			model.addConstr(rcapjs(j, s) == (bcj[j] - sum(j + 1, J, [&](int i) { return njs(i, s) * cj(i); })));

			GRBVar kjsarr[] = { kjs(j,s) };
			model.addGenConstrMin(njs(j, s), kjsarr, 1, scenarios(s, j));
		}
	}

	return {bcj, njs};
}

BookingAcceptVars modelBuilderForEconomiesOfScale(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) {
	int J = sim.getNumClasses(),
		C = sim.getC(),
		S = scenarios.getM(),
		U = C + 1;

	auto rj = [&](int j) { return sim.getCustomer(j).revenuePerReq; };

	auto cju = [&](int j, int u) {
		return ((MultiClassSimulation *)&sim)->eosConsumption(j, u);
	};

	vector<GRBVar> bcj = Helpers::constructVector<GRBVar>(J, [&](int j) { return model.addVar(j == 0 ? C : 0, C, 0.0, GRB_INTEGER, "bcj" + to_string(j)); });

	vector<Matrix<GRBVar>> njsu = Helpers::constructVector<Matrix<GRBVar>>(J, [&](int j) {
		Matrix<GRBVar> nsu(S, U, [&](int s, int u) {
			string caption = "njsu" + to_string(j) + "," + to_string(s) + "," + to_string(u);
			return model.addVar(0.0, 1.0, 0.0, GRB_BINARY, caption);
		});
		return nsu;
	});

	vector<Matrix<GRBVar>> nbjsu = Helpers::constructVector<Matrix<GRBVar>>(J, [&](int j) {
		Matrix<GRBVar> nbsu(S, U, [&](int s, int u) {
			string caption = "nbjsu" + to_string(j) + "," + to_string(s) + "," + to_string(u);
			return model.addVar(0.0, 1.0, 0.0, GRB_BINARY, caption);
		});
		return nbsu;
	});

	Matrix<GRBVar> njs(J, S, [&](int j, int s) { return model.addVar(0, C, 0.0, GRB_INTEGER, "njs" + to_string(j) + "," + to_string(s)); });
	Matrix<GRBVar> nbjs(J, S, [&](int j, int s) { return model.addVar(0, C, 0.0, GRB_INTEGER, "nbjs" + to_string(j) + "," + to_string(s)); });

	//model.setObjective(1.0 / (double)S * sum3D(J, S ,U, [&](int j, int s, int u) {  return njsu[j](s, u) * u * rj(j); }), GRB_MAXIMIZE);
	model.setObjective(1.0 / (double)S * sum2D(J, S, [&](int j, int s) {  return njs(j, s) * rj(j); }), GRB_MAXIMIZE);

	for (int j = 0; j < J; j++) {
		if (j + 1 < J)
			model.addConstr(bcj[j] >= bcj[j + 1]);

		for (int s = 0; s < S; s++) {
			model.addConstr(sum(U, [&](int u) { return njsu[j](s, u); }) == 1);
			model.addConstr(sum(U, [&](int u) { return nbjsu[j](s, u); }) == 1);

			model.addConstr(sum(U, [&](int u) { return njsu[j](s, u) * u; }) == njs(j, s));
			model.addConstr(sum(U, [&](int u) { return nbjsu[j](s, u) * u; }) == nbjs(j, s));

			GRBVar nbjsarr[] = { nbjs(j,s) };
			model.addGenConstrMin(njs(j, s), nbjsarr, 1, scenarios(s, j));

			model.addConstr(sum2D(j + 1, J, 0, U, [&](int i, int u) { return njsu[i](s, u) * u * cju(i, u); }) + sum(U, [&](int u) { return nbjsu[j](s, u) * u * cju(j, u); }) <= bcj[j]);
			model.addConstr(sum2D(j + 1, J, 0, U, [&](int i, int u) { return njsu[i](s, u) * u * cju(i, u); }) + sum(U, [&](int u) { return nbjsu[j](s, u) * (u + 1) * cju(j, u + 1); }) >= bcj[j] + globals::EPSILON2);
		}
	}

	//cout << "EPISLON2: " << globals::EPSILON2 << endl;

	return {bcj, njs};
}

BookingAcceptVars modelBuilderForNewFormulation(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) {
	int J = sim.getNumClasses(),
			C = sim.getC(),
			S = scenarios.getM();

	auto cj = [&](int j) { return sim.getCustomer(j).consumptionPerReqMean; };
	auto rj = [&](int j) { return sim.getCustomer(j).revenuePerReq; };

	vector<GRBVar> bcj = Helpers::constructVector<GRBVar>(J, [&](int j) { return model.addVar(j == 0 ? C : 0, C, 0.0, GRB_INTEGER, "bcj" + to_string(j)); });

	Matrix<GRBVar> njs(J, S, [&](int j, int s) {
		string caption = "njs" + to_string(j) + "," + to_string(s);
		return model.addVar(0.0, C, 0.0, GRB_INTEGER, caption);
	});

	Matrix<GRBVar> nbjs(J, S, [&](int j, int s) {
		string caption = "nbjs" + to_string(j) + "," + to_string(s);
		return model.addVar(0.0, C, 0.0, GRB_INTEGER, caption);
	});

	model.setObjective(1.0 / (double)S * sum2D(J, S, [&](int j, int s) {  return njs(j, s) * rj(j); }), GRB_MAXIMIZE);

	for (int j = 0; j < J; j++) {
		if (j + 1 < J)
			model.addConstr(bcj[j] >= bcj[j + 1]);

		for (int s = 0; s < S; s++) {
			GRBVar nbjsarr[] = { nbjs(j,s) };
			model.addGenConstrMin(njs(j,s), nbjsarr, 1, scenarios(s,j));

			model.addConstr(sum(j+1, J, [&](int i) { return njs(i,s) * cj(i); }) + nbjs(j,s) * cj(j) <= bcj[j]);
			model.addConstr(sum(j+1, J, [&](int i) { return njs(i,s) * cj(i); }) + (nbjs(j,s) + 1.0) * cj(j) >= bcj[j] + globals::EPSILON);
		}
	}

	return {bcj, njs};
}

void addConditionalValueAtRiskConsiderationToModel(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const Matrix<GRBVar> &njs) {
	int J = sim.getNumClasses(),
		S = scenarios.getM();

	const double cvarUpperBound = ceil(sim.getC() / sim.getCustomer(0).consumptionPerReqMean) * sim.getCustomer(0).revenuePerReq;
	const GRBVar cvar = model.addVar(0.0, cvarUpperBound, 0.0, GRB_CONTINUOUS, "CVaR");

	const GRBVar omega0 = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "omega0");
	vector<GRBVar> omega_s = Helpers::constructVector<GRBVar>(S, [&](int s) { return model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "omega_s=" + to_string(s)); });

	const double profitWeight = sim.getPsi();
	const double alpha = sim.getAlpha();

	vector<GRBVar> G_s = Helpers::constructVector<GRBVar>(S, [&](int s) { return model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "G_s=" + to_string(s));  });

	model.setObjective(profitWeight * 1.0 / (double)S * sum(S, [&](int s) {  return G_s[s]; }) + (1 - profitWeight) * cvar, GRB_MAXIMIZE);
	model.addConstr(cvar == omega0 - 1 / (1 - alpha) * sum(S, [&omega_s](int s) { return omega_s[s]; }));
	for (int s = 0; s < S; s++) {
		model.addConstr(G_s[s] == sum(J, [&](int j) { return sim.getCustomer(j).revenuePerReq * njs(j, s); }));
		model.addConstr(omega0 - omega_s[s] <= G_s[s]);
	}
}

BookingAcceptVars modelBuilderForConditionalValueAtRisk(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) {
	auto pair = modelBuilderForNewFormulation(sim, model, scenarios, {});
	addConditionalValueAtRiskConsiderationToModel(sim, model, scenarios, pair.njs);
	return pair;
}

BookingAcceptVars modelBuilderForStochasticConsumptions(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) {
	if (!consumptionScenarioFunc.is_initialized()) {
		throw runtime_error("No consumption scenarios provided for stochastic consumption model!");
	}

	int J = sim.getNumClasses(),
		C = sim.getC(),
		S = scenarios.getM();

	auto cjs = [&](int j, int s) { return (consumptionScenarioFunc.get())(0, s, j); };
	auto rj = [&](int j) { return sim.getCustomer(j).revenuePerReq; };

	vector<GRBVar> bcj = Helpers::constructVector<GRBVar>(J, [&](int j) { return model.addVar(j == 0 ? C : 0, C, 0.0, GRB_INTEGER, "bcj" + to_string(j)); });

	Matrix<GRBVar> njs(J, S, [&](int j, int s) {
		string caption = "njs" + to_string(j) + "," + to_string(s);
		return model.addVar(0.0, C, 0.0, GRB_INTEGER, caption);
	});

	Matrix<GRBVar> nbjs(J, S, [&](int j, int s) {
		string caption = "nbjs" + to_string(j) + "," + to_string(s);
		return model.addVar(0.0, C, 0.0, GRB_INTEGER, caption);
	});

	model.setObjective(1.0 / (double)S * sum2D(J, S, [&](int j, int s) {  return njs(j, s) * rj(j); }), GRB_MAXIMIZE);

	for (int j = 0; j < J; j++) {
		if (j + 1 < J)
			model.addConstr(bcj[j] >= bcj[j + 1]);

		for (int s = 0; s < S; s++) {
			GRBVar nbjsarr[] = { nbjs(j,s) };
			model.addGenConstrMin(njs(j, s), nbjsarr, 1, scenarios(s, j));

			model.addConstr(sum(j + 1, J, [&](int i) { return njs(i, s) * cjs(i, s); }) + nbjs(j, s) * cjs(j, s) <= bcj[j]);
			model.addConstr(sum(j + 1, J, [&](int i) { return njs(i, s) * cjs(i, s); }) + (nbjs(j, s) + 1.0) * cjs(j, s) >= bcj[j] + globals::EPSILON);
		}
	}

	return {bcj, njs};
}

BookingAcceptVars modelBuilderForStochasticEconomiesOfScale(const AbstractSimulation &sim, GRBModel &model, const DemandScenarioList &scenarios, const boost::optional<ConsumptionScenarioFunc&> consumptionScenarioFunc) {
	int J = sim.getNumClasses(),
			C = sim.getC(),
			S = scenarios.getM(),
			U = C + 1;

	auto rj = [&](int j) { return sim.getCustomer(j).revenuePerReq; };

    auto cjsu = [&](int j, int s, int u) { return (consumptionScenarioFunc.get())(u, s, j); };

	vector<GRBVar> bcj = Helpers::constructVector<GRBVar>(J, [&](int j) { return model.addVar(j == 0 ? C : 0, C, 0.0, GRB_INTEGER, "bcj" + to_string(j)); });

	vector<Matrix<GRBVar>> njsu = Helpers::constructVector<Matrix<GRBVar>>(J, [&](int j) {
		Matrix<GRBVar> nsu(S, U, [&](int s, int u) {
			string caption = "njsu" + to_string(j) + "," + to_string(s) + "," + to_string(u);
			return model.addVar(0.0, 1.0, 0.0, GRB_BINARY, caption);
		});
		return nsu;
	});

	vector<Matrix<GRBVar>> nbjsu = Helpers::constructVector<Matrix<GRBVar>>(J, [&](int j) {
		Matrix<GRBVar> nbsu(S, U, [&](int s, int u) {
			string caption = "nbjsu" + to_string(j) + "," + to_string(s) + "," + to_string(u);
			return model.addVar(0.0, 1.0, 0.0, GRB_BINARY, caption);
		});
		return nbsu;
	});

	Matrix<GRBVar> njs(J, S, [&](int j, int s) { return model.addVar(0, C, 0.0, GRB_INTEGER, "njs" + to_string(j) + "," + to_string(s)); });
	Matrix<GRBVar> nbjs(J, S, [&](int j, int s) { return model.addVar(0, C, 0.0, GRB_INTEGER, "nbjs" + to_string(j) + "," + to_string(s)); });

	//model.setObjective(1.0 / (double)S * sum3D(J, S ,U, [&](int j, int s, int u) {  return njsu[j](s, u) * u * rj(j); }), GRB_MAXIMIZE);
	model.setObjective(1.0 / (double)S * sum2D(J, S, [&](int j, int s) {  return njs(j, s) * rj(j); }), GRB_MAXIMIZE);

	for (int j = 0; j < J; j++) {
		if (j + 1 < J)
			model.addConstr(bcj[j] >= bcj[j + 1]);

		for (int s = 0; s < S; s++) {
			model.addConstr(sum(U, [&](int u) { return njsu[j](s, u); }) == 1);
			model.addConstr(sum(U, [&](int u) { return nbjsu[j](s, u); }) == 1);

			model.addConstr(sum(U, [&](int u) { return njsu[j](s, u) * u; }) == njs(j, s));
			model.addConstr(sum(U, [&](int u) { return nbjsu[j](s, u) * u; }) == nbjs(j, s));

			GRBVar nbjsarr[] = { nbjs(j,s) };
			model.addGenConstrMin(njs(j, s), nbjsarr, 1, scenarios(s, j));

			model.addConstr(sum2D(j + 1, J, 0, U, [&](int i, int u) { return njsu[i](s, u) * u * cjsu(i, s, u); }) + sum(U, [&](int u) { return nbjsu[j](s, u) * u * cjsu(j, s, u); }) <= bcj[j]);
			model.addConstr(sum2D(j + 1, J, 0, U, [&](int i, int u) { return njsu[i](s, u) * u * cjsu(i, s, u); }) + sum(U, [&](int u) { return nbjsu[j](s, u) * (u + 1) * cjsu(j, s, u + 1); }) >= bcj[j] + globals::EPSILON2);
		}
	}

	//cout << "EPISLON2: " << globals::EPSILON2 << endl;

	return { bcj, njs};
}