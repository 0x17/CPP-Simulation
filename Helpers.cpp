//
// Created by André Schnabel on 30.10.16.
//

#include "Helpers.h"
#include "Stopwatch.h"

#include <boost/math/distributions/normal.hpp>
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

static mt19937 *gen = new mt19937(42);

void Helpers::resetSeed(int seed) {
	if(gen) {
		delete gen;
	}
	gen = new mt19937(seed);
}

double Helpers::pickNormal(double mean, double stddev) {
	normal_distribution<double> d(mean, stddev);
    return d(*gen);
}

double Helpers::invNormal(double x, double mean, double stddev) {
	boost::math::normal d(mean, stddev);
	return quantile(d, x);
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

namespace Helpers {
	Tracer::Tracer(const string &filePrefix) : f(filePrefix + ".txt") {
		sw.start();
		lupdate = chrono::system_clock::now();
		last_slvtime = 0.0;
		if(!f.is_open())
			throw runtime_error("Unable to create " + filePrefix + ".txt!");
		f << "slvtime;bks_objval\n";
		trace(0.0, 0.0f);
	}

	Tracer::~Tracer() {
		f.close();
	}

	void Tracer::trace(double slvtime, double bks_objval, bool trunc_secs) {
		double insecs = (slvtime / 1000.0);
		if (trunc_secs) insecs = trunc(insecs);
		f << (boost::format("%.2f") % insecs) << ";" << bks_objval << endl;
	}

	void Tracer::intervalTrace(double bks_objval) {
		double slvtime = sw.look();
		double deltat = chrono::duration<double, milli>(chrono::system_clock::now() - lupdate).count();
		if(slvtime < 1000.0 && deltat >= MSECS_BETWEEN_TRACES_SHORT) {
			lupdate = chrono::system_clock::now();
			trace(slvtime, bks_objval);
		} else if(slvtime >= 1000.0 && last_slvtime < 1000.0) {
			lupdate = chrono::system_clock::now();
			trace(slvtime, bks_objval, true);
		} else if(slvtime >= 1000.0 && deltat >= MSECS_BETWEEN_TRACES_LONG) {
			//cout << "Nodes visited = " << nodeCtr << ", Boundings = " << boundCtr << ", Opt = " << lb << ", Time = " << (boost::format("%.2f") % (sw.look() / 1000.0)) << endl;
			lupdate = chrono::system_clock::now();
			trace(slvtime, bks_objval, true);
		}
		last_slvtime = slvtime;
	}
}