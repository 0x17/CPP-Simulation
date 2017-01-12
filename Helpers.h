//
// Created by André Schnabel on 30.10.16.
//

#ifndef CPP_SIMULATION_HELPERS_H
#define CPP_SIMULATION_HELPERS_H

#include <string>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <fstream>
#include "Stopwatch.h"

namespace Helpers {
    std::string slurp(const std::string& filename);

	void resetSeed(int seed);
    double pickNormal(double mean, double stddev);
	double invNormal(double x, double mean, double stddev);
	double pickNormalDescriptive(double mean, double stddev, int scenarioIx, int nscenarios);

	double vecAverage(const std::vector<double> & nums);

	template<class A, class Func>
	std::vector<A> constructVector(int size, Func f) {
		std::vector<A> v(size);
		for (int i = 0; i<size; i++) {
			v[i] = f(i);
		}
		return v;
	}

	void spit(const std::string &s, const std::string &filename);
	void spitAppend(const std::string &s, const std::string &filename);

	inline int randRangeIncl(int lb, int ub) {
		return lb + rand() % (ub - lb + 1);
	}

	inline float randUnitFloat() {
		return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
	}

	inline double randUnitDouble() {
		return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
	}

	class Tracer {
		std::ofstream f;
		std::chrono::time_point<std::chrono::system_clock> lupdate;
		double last_slvtime;
		Stopwatch sw;
	public:
		Tracer(const std::string &filePrefix = "SolverTrace");
		~Tracer();
		void trace(double slvtime, float bks_objval, bool trunc_secs = false);
		void intervalTrace(float bks_objval);
	};
};


#endif //CPP_SIMULATION_HELPERS_H
