//
// Created by André Schnabel on 30.10.16.
//

#ifndef CPP_SIMULATION_HELPERS_H
#define CPP_SIMULATION_HELPERS_H

#include <string>
#include <vector>

namespace Helpers {
    std::string slurp(const std::string& filename);

	void resetSeed(int seed);
    double pickNormal(double mean, double stddev);
	double invNormal(double x, double mu, double sigma);

	double vecAverage(const std::vector<double> & nums);

	template<class A, class Func>
	std::vector<A> constructVector(int size, Func f) {
		std::vector<A> v(size);
		for (int i = 0; i<size; i++) {
			v[i] = f(i);
		}
		return v;
	}
};


#endif //CPP_SIMULATION_HELPERS_H
