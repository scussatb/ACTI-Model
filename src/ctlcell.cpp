#include "ctlcell.h"
#include "tumorcell.h"

#include <limits>

double CTLCell::flySpeed = 8.66;
std::exponential_distribution<double> CTLCell::flyExpDist(1.0);
std::normal_distribution<double> CTLCell::flyBrownianDistanceNormDist(flySpeed, 2.5);
std::uniform_real_distribution<double> CTLCell::flyAngleUniDist(0., 2.0 * M_PI);
pareto_distribution CTLCell::flyLengthParDist(1.99891, 1);

double CTLCell::inhibitoryProbability = 0.00275; // per minute

double CTLCell::hitStrength = 1.0 / 55.0;

CTLCell::CTLCell(MecaCell::Vec v):	
	VCell(v) {
	type = CTL;
	age = 0.;
	nextChange = 0.;
	nChange=0;

	inhibitorySignal = 0.0;
	
	state = flying;

	color[0] = 0.0;
	color[1] = 0.7;
	color[2] = 0.0;

	nbKilledCell = 0;
	tumorSynapticTime = 0;
	ctlSynapticTime = 0;
}

CTLCell::~CTLCell() {
}

double CTLCell::getAdhesionWith(const VCell *c) {
	return (c->type == Tumor)?0.8:0.05;
}

void CTLCell::randomMove(double dt) {
	// check if the cell is outside the box
	if (position.x * position.x + position.z * position.z > 0.9*0.9*Scenario::BOX_SIZE * Scenario::BOX_SIZE/4.) {
		double speed = flyDirection.length();
		MecaCell::Vec directionToCenter(-position.x, 0.0, -position.z);
		directionToCenter.normalize();
		directionToCenter.y = 0;
	    setVelocity(MecaCell::Vec(directionToCenter * speed));
		return;
	}
	color[0] = 0.0;
	color[1] = 0.7 * (1.0 - inhibitorySignal);
	color[2] = 0.7 * inhibitorySignal;
	if (age > nextChange) {
		double angle = flyAngleUniDist(MecaCell::globalRand);
		flyDirection.x = std::cos(angle);
		flyDirection.y = 0.0;
		flyDirection.z = std::sin(angle);
		flyDirection.normalize();
		flyDirection.y = 0.0;
		flyDirection *= flyBrownianDistanceNormDist(MecaCell::globalRand) * Scenario::MinutesPerTick * dt;

		nextChange = age + Scenario::MinutesPerTick * dt ;  
		nChange++;
	}
	setVelocity(flyDirection);
}

VCell *CTLCell::fly(double dt) {
	randomMove(dt);
	return nullptr;
}

VCell *CTLCell::touchingTumor(double dt) {
	int nTumorInContact = 0;
	for (auto *c : getConnectedCells()) {
		if (c->type == Tumor) {
			nTumorInContact++;
			TumorCell* tc = (TumorCell*)c;
			if (tc->isDetected()) {
				color[0] = 0.0;
				color[1] = 0.7 * (1.0 - inhibitorySignal);
				color[2] = 0.7 * inhibitorySignal;
				// a tumor cell has been detected => attacking
				if (std::find(scenario->collidingCTLCells.begin(), scenario->collidingCTLCells.end(), this) == scenario->collidingCTLCells.end()) {
					// adding myself as attractor
					scenario->collidingCTLCells.push_back(this);
				}
			} else {
				color[0] = 0.0;
				color[1] = 0.7 * (1.0 - inhibitorySignal);
				color[2] = 0.7 * inhibitorySignal;
			}
		}
	}

	// no more connection => get back to flying
	if (nTumorInContact == 0) {
		state = flying;
		return nullptr;
	}

	return nullptr;	
}

VCell *CTLCell::updateBehavior(double dt) {
	if (isDead()) {
		return nullptr;
	}
	if (getVelocity().y<-10.) {
		// currently falling, do nothing
		return nullptr;
	}
	age+=dt;
	receiveForce(0.01, MecaCell::Vec(-position.x, 0, -position.z), false);
	if (Scenario::synapticContacts) {
	   	for (auto *c : getConnectedCells()) {
			if (c->type == Tumor) tumorSynapticTime++;
			else  if (c->type == CTL) ctlSynapticTime++;
		}
	}
	
	switch (state) {
		case flying:
			return fly(dt);
			break;
		case contactWithTumor:
			return touchingTumor(dt);
			break;
	}
	return nullptr;
}


