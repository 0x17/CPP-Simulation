#include <boost/algorithm/clamp.hpp>
#include <functional>
#include <cmath>
#include <iostream>

#include "PSSolver.h"
#include "Matrix.h"
#include "Helpers.h"
#include "Stopwatch.h"
#include "Globals.h"

using namespace std;

class Swarm {
public:
	Swarm(int _swarmSize, int _numClasses, int _C, function<double(vector<int>)> _objective, boost::optional<vector<int>> seedSolution);	
	void update();
	Result getBestResult() const;

	void writeSwarmToFile(const string &filename) const;

private:
	void initStartingPositions();	
	void initLocalGlobalBests();
	void initStartingVelocities();

	void restrictParticleToBounds(int particleIndex);
	void updateParticlePositionAndVelocity(int particleIndex);
	void updateLocalGlobalBests(int particleIndex);

	int swarmSize, numClasses, C;
	Matrix<int> particles, personalBests;
	Matrix<double> velocities;
	vector<double> personalBestObjectives;
	double globalBestObjective;
	vector<int> globalBest;
	function<double(vector<int>)> objective;
	boost::optional<vector<int>> seedSolution;
};

PSSolver::PSSolver(const AbstractSimulation &_sim) : BookingLimitOptimizer("ParticleSwarm", _sim) {}

Result PSSolver::solve(vector<vector<int>>& scenarios) {
	const int	iterlimit = -1,
				swarmSize = 20;
	const double timelimit = globals::TIME_LIMIT;

	auto objective = [&](vector<int> bookingLimits) {
		return sim.averageRevenueOfSimulation(bookingLimits, scenarios);
	};

	Swarm s(swarmSize, sim.getNumClasses(), sim.getC(), objective, heuristicBookingLimits);
	Stopwatch sw;

	sw.start();
	double tstart = sw.look();
	for (int i = 0; (iterlimit != -1 && i < iterlimit) || (timelimit != -1.0 && sw.look() - tstart < timelimit * 1000.0); i++) {
		s.update();
#ifndef __APPLE__
		cout << "Particle swarm iteration " << i << "\r" << flush;
#endif
		//s.writeSwarmToFile("swarmIteration" + to_string(i + 1) + ".txt");
	}

	return s.getBestResult();
}

Swarm::Swarm(int _swarmSize, int _numClasses, int _C, function<double(vector<int>)> _objective, boost::optional<vector<int>> _seedSolution) :
	swarmSize(_swarmSize),
	numClasses(_numClasses),
	C(_C),
	particles(swarmSize, _numClasses),
	personalBests(swarmSize, _numClasses),
	velocities(swarmSize, _numClasses),
	personalBestObjectives(_swarmSize),
	globalBestObjective(numeric_limits<double>::lowest()),
	globalBest(_numClasses),
	objective(_objective),
	seedSolution(_seedSolution) {

	initStartingPositions();
	initLocalGlobalBests();
	initStartingVelocities();
}

void Swarm::restrictParticleToBounds(int particleIndex) {
	particles(particleIndex, 0) = boost::algorithm::clamp(particles(particleIndex, 0), 0, C);
	for(int j=1; j<numClasses; j++) {
		particles(particleIndex, j) = boost::algorithm::clamp(particles(particleIndex, j), 0, particles(particleIndex, j - 1));
	}
}

void Swarm::updateParticlePositionAndVelocity(int particleIndex) {
	const double omega = 0.5, phi_p = 0.2, phi_g = 0.4;

	double	r_p = Helpers::randUnitDouble(),
			r_g = Helpers::randUnitDouble();

	for(int j=1; j<numClasses; j++) {			
		velocities(particleIndex, j) = omega * velocities(particleIndex, j) + phi_p * r_p * (personalBests(particleIndex, j) - particles(particleIndex, j)) + phi_g * r_g * (globalBest[j] - particles(particleIndex, j));
		particles(particleIndex, j) = (int)floor(particles(particleIndex, j) + velocities(particleIndex, j));
	}
}

void Swarm::updateLocalGlobalBests(int particleIndex) {
	double newObj = objective(particles.row(particleIndex));
	if(newObj > personalBestObjectives[particleIndex]) {
		// update personal best
		personalBestObjectives[particleIndex] = newObj;
		for(int j=1; j<numClasses; j++) {
			personalBests(particleIndex, j) = particles(particleIndex, j);
		}			
		// update global best
		if(newObj > globalBestObjective) {
			globalBestObjective = newObj;
			globalBest = particles.row(particleIndex);
		}
	}
}

void Swarm::update() {
	static Helpers::Tracer tr("ParticleSwarmTrace");

	for (int i = 0; i < swarmSize; i++) {
		updateParticlePositionAndVelocity(i);
		restrictParticleToBounds(i);
		updateLocalGlobalBests(i);
	}

	tr.intervalTrace(globalBestObjective);
}

Result Swarm::getBestResult() const {
	return{ globalBest, globalBestObjective };
}

void Swarm::initStartingPositions() {
	for (int i = 0; i < swarmSize; i++) {
		particles(i, 0) = C;
		for (int j = 1; j < numClasses; j++) {
			particles(i, j) = Helpers::randRangeIncl(0, particles(i, j-1));
		}
	}

	// ensure heuristic policy is included in swarm
	if (seedSolution) {
		for (int j = 0; j < numClasses; j++)
			particles(0, j) = (*seedSolution)[j];
	}
}

void Swarm::writeSwarmToFile(const string &filename) const {
	if(numClasses != 3) return;

	Helpers::spit("b2;b3;obj\n", filename);
	for(int i=0; i<swarmSize; i++) {
		Helpers::spitAppend(to_string(particles(i,1)) + ";" + to_string(particles(i,2)) + ";" + to_string(objective(particles.row(i))) + "\n", filename);
	}
}

void Swarm::initLocalGlobalBests() {
	personalBests = particles;
	for (int i = 0; i < swarmSize; i++) {
		double obj = objective(personalBests.row(i));
		personalBestObjectives[i] = obj;
		if (obj > globalBestObjective) {
			globalBest = personalBests.row(i);
			globalBestObjective = obj;
		}
	}

	//writeSwarmToFile("initialSwarm.txt");
}

void Swarm::initStartingVelocities() {
	velocities.foreachAssign([](int i, int j) { return 0; });
}
