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
#include <list>
#include <memory>
#include "Stopwatch.h"
#include "json11.hpp"

namespace Helpers {

	template<class T>
	void swap(std::vector<T> &v, int ix1, int ix2) {
		T tmp = v[ix1];
		v[ix1] = v[ix2];
		v[ix2] = tmp;
	}

    std::string slurp(const std::string& filename);
	json11::Json readJsonFromFile(const std::string &filename);

	void resetSeed(int seed);
    double pickNormal(double mean, double stddev);
	double invNormal(double x, double mean, double stddev);
	double pickNormalDescriptive(double mean, double stddev, int scenarioIx, int nscenarios);

	std::vector<double> generateNormalDistributionDescriptiveSamplingLUT(int sampleSize, double mean, double stddev);

	double pickNextWithLUT(std::vector<double> &lut, int& drawnCounter);
	double pickNextWithLUT(std::vector<double> &lut);

	double vecAverage(const std::vector<int> & nums);
	double vecAverage(const std::vector<double> & nums);
	double vecAverageSubRange(const std::vector<double> &nums, int startIndexInclusively, int endIndexExclusively);

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
		std::unique_ptr<std::ofstream> f;
		std::chrono::time_point<std::chrono::system_clock> lupdate;
		double last_slvtime;
		Stopwatch sw;
	public:
		explicit Tracer(const std::string &filePrefix = "SolverTrace");
		~Tracer();
		void trace(double slvtime, double bks_objval, bool trunc_secs = false);
		void intervalTrace(double bks_objval);
	};

	std::list<std::string> extractArguments(int argc, const char **argv);
};


#endif //CPP_SIMULATION_HELPERS_H
