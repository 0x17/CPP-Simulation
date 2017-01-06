#include "PSSolver.h"
#include "Matrix.h"
#include "Helpers.h"
#include <functional>

class Swarm {
public:
	Swarm(int _swarmSize, int _numClasses, int _C, std::function<double(std::vector<int>)> _objective);
	void update();
	Result getBestResult() const;

private:
	void initStartingPositions();
	void initLocalGlobalBests();
	void initStartingVelocities();

	int swarmSize, numClasses, C;
	Matrix<int> particles, personalBests, velocities;
	double globalBestObjective;
	std::vector<int> globalBest;
	std::function<double(std::vector<int>)> objective;
};

PSSolver::PSSolver(AbstractSimulation &_sim) : BookingLimitOptimizer("ParticleSwarm", _sim) {}

Result PSSolver::solve(std::vector<std::vector<int>>& scenarios) {
	const int iterlimit = 10;
	const int swarmSize = 100;

	auto objective = [&](std::vector<int> bookingLimits) {
		return Helpers::vecAverage(sim.runSimulation(bookingLimits, scenarios));
	};

	Swarm s(swarmSize, sim.getNumClasses(), sim.getC(), objective);
	
	for (int i = 0; i < iterlimit; i++)
		s.update();

	return s.getBestResult();
}

Swarm::Swarm(int _swarmSize, int _numClasses, int _C, std::function<double(std::vector<int>)> _objective) :
	swarmSize(_swarmSize),
	numClasses(_numClasses),
	C(_C),
	particles(swarmSize, _numClasses),
	personalBests(swarmSize, _numClasses),
	velocities(swarmSize, _numClasses),
	objective(_objective) {

	initStartingPositions();
	initLocalGlobalBests();
	initStartingVelocities();
}

void Swarm::update() {
	const float omega = 0.2, phi_p = 0.2, phi_g = 0.2;

	for (int i = 0; i < swarmSize; i++) {
		int r_p = Helpers::randUnitFloat(),
			r_g = Helpers::randUnitFloat();

		for(int j=0; j<numClasses; j++) {			
			velocities(i, j) = omega * velocities(i, j) + phi_p * r_p * (personalBests(i, j) - particles(i, j)) + phi_g * r_g * (globalBest[j] - particles(i, j));
			particles(i, j) = particles(i, j) + velocities(i, j);
		}

		if(objective(particles.row(i)) > objective(personalBests.row(i))) {
			for(int j=0; j<numClasses; j++) {
				particles(i, j) = personalBests(i, j);
			}			
			if(objective(particles.row(i)) > globalBestObjective) {
				globalBestObjective = objective(particles.row(i));
				globalBest = particles.row(i);
			}
		}
	}
}

Result Swarm::getBestResult() const {
	return{ globalBest, globalBestObjective };
}

void Swarm::initStartingPositions() {
	for (int i = 0; i < swarmSize; i++) {
		particles(i, 0) = Helpers::randRangeIncl(0, C);
		for (int j = 1; j < numClasses; j++) {
			particles(i, j) = Helpers::randRangeIncl(0, particles(i, j));
		}
	}
}

void Swarm::initLocalGlobalBests() {
	personalBests = particles;
	std::vector<int> globalBest(numClasses);
	globalBestObjective = std::numeric_limits<int>::min();
	for (int i = 0; i < swarmSize; i++) {
		double obj = objective(personalBests.row(i));
		if (obj > globalBestObjective) {
			globalBest = personalBests.row(i);
			globalBestObjective = obj;
		}
	}
}

void Swarm::initStartingVelocities() {
	velocities.foreachAssign([](int i, int j) { return 0; });
}