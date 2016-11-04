//
// Created by André Schnabel on 30.10.16.
//

#include "Helpers.h"

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

double Helpers::invNormal(double x, double mu, double sigma) {
	return 0.5 * erfc((x - mu) / (sigma * sqrt(2)));
}

double Helpers::vecAverage(const vector<double>& nums) {
	return accumulate(nums.begin(), nums.end(), 0.0) / (double)nums.size();
}
