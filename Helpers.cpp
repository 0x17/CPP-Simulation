//
// Created by André Schnabel on 30.10.16.
//

#include "Helpers.h"

#include <boost/math/special_functions/erf.hpp>
#include <fstream>
#include <random>
#include <map>
#include <numeric>
using namespace std;

string Helpers::slurp(const string &filename) {
    ifstream ifs(filename);
    string content((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
    return content;
}

static mt19937 *gen = nullptr;

void Helpers::resetSeed(int seed) {
	if(gen) {
		delete gen;
	}
	gen = new mt19937(seed);
}

double Helpers::pickNormal(double mean, double stddev) {
    static map<pair<double, double>, normal_distribution<double>> distCache;
    auto params = make_pair(mean, stddev);
    auto it = distCache.find(params);
    if(it == distCache.end()) {
        normal_distribution<double> d(mean, stddev);
        distCache[params] = d;
    }
    return distCache[params](*gen);
}

double Helpers::invNormal(double x, double mean, double stddev) {
	return mean + sqrt(2) * stddev * boost::math::erfc_inv(2.0 * x);
}

double Helpers::pickNormalDescriptive(double mean, double stddev, int scenarioIx, int nscenarios) {
	return invNormal(((double)scenarioIx + 0.5) / (double)nscenarios, mean, stddev);
}

double Helpers::vecAverage(const vector<double>& nums) {
	return accumulate(nums.begin(), nums.end(), 0.0) / (double)nums.size();
}

void Helpers::spit(const std::string &s, const std::string &filename) {
	std::ofstream f(filename);
	if(f.is_open()) {
		f << s;
		f.close();
	}
}

void Helpers::spitAppend(const string &s, const string &filename) {
	ofstream f(filename, ios_base::app);
	if (f.is_open()) {
		f << s;
		f.close();
	}
}