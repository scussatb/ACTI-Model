#include "ctlcell.h"
#include "tumorcell.h"

#include <limits>

double CTLCell::attractionRadius = 0.0;

double CTLCell::attractionMoveSpeed = 8.66;
std::exponential_distribution<double> CTLCell::flyExpDist(1.0);
std::normal_distribution<double> CTLCell::flyBrownianDistanceNormDist(attractionMoveSpeed, 2.5);
std::uniform_real_distribution<double> CTLCell::flyAngleUniDist(0., 2.0 * M_PI);
pareto_distribution CTLCell::flyLengthParDist(1.99891, 1);
std::exponential_distribution<double> CTLCell::flyJumpDurationExpDist(4.927135);
double CTLCell::tScanMax = 5;
std::normal_distribution<double> CTLCell::hitStrengthNormDist(1.0/55.0, 0.0); // per minute
double CTLCell::inhibitoryRecoverySpeed = 0.00275; // per minute
double CTLCell::attractingReduceSpeed = 1.0; // per minute

CTLCell::CTLCell(MecaCell::Vec v):	
	VCell(v) {
	type = CTL;
	age = 0.;
	nextChange = 0.;
	nChange=0;

	inhibitorySignal = 0.0;
	hitStrength = std::max(0.0, hitStrengthNormDist(MecaCell::globalRand));
	
	state = flying;

	color[0] = 0.0;
	color[1] = 0.7;
	color[2] = 0.0;

	nbKilledCell = 0;
	tumorSynapticTime = 0;
	ctlSynapticTime = 0;
	attractingCounter = 0.0;
}

CTLCell::~CTLCell() {
}

double CTLCell::getAdhesionWith(const VCell *c) {
	return (c->type == Tumor)?0.8:0.05;
}

void CTLCell::randomMove(double dt) {
	static double lastX = 0;
	static double lastY = 0;
	static double lastZ = 0;

	// check if the cell is outside the box
	if (position.x * position.x + position.z * position.z > 0.9*0.9*Scenario::BOX_SIZE * Scenario::BOX_SIZE/4.) {
		double speed = flyDirection.length();
		MecaCell::Vec directionToCenter(-position.x, 0.0, -position.z);
		directionToCenter.normalize();
		directionToCenter.y = 0;
	        setVelocity(MecaCell::Vec(directionToCenter * speed));
		return;
	}
	if (flyType == brownian) {
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
			MecaCell::Vec centerAttrac(-position.x, 0, -position.z);
			centerAttrac.normalize();
			centerAttrac *= 1. * Scenario::MinutesPerTick * dt;
			flyDirection += centerAttrac;

			nextChange = age + Scenario::MinutesPerTick * dt ;  
			nChange++;
		}
		setVelocity(flyDirection);
	} else { // jump
		color[0] = 0.0;
		color[1] = 0.7 * (1.0 - inhibitorySignal);
		color[2] = 0.7 * inhibitorySignal;
		setVelocity(flyDirection);
		if (age>flyTimeMode) {
			// end of the jump
			flyType = brownian;
			flyTimeMode += flyJumpDurationExpDist(MecaCell::globalRand);
		}
	}
}

VCell *CTLCell::fly(double dt) {
	// first move
	randomMove(dt);

	// looking for connected tumor cells
	if (tScan < tScanMax) {
		for (auto *c : getConnectedCells()) {
			if (c->type == Tumor) {
				state = contactWithTumor;
				setVelocity(MecaCell::Vec(0,0,0));
				tScan = 0;
				return nullptr;
			}
		}
	} else {
		tScan -= dt * Scenario::MinutesPerTick;
	}

	// looking for attractions
	attractionDirection.x = 0;
	attractionDirection.y = 0;
	attractionDirection.z = 0;
	double minDist = attractionRadius * getRadius() * attractionRadius * getRadius();
	CTLCell *closerCTL = nullptr;
	for (auto *c : scenario->collidingCTLCells) {
		if (c != this) {
			double dist = (position.x - c->position.x) * (position.x - c->position.x) + 
				(position.y - c->position.y) * (position.y - c->position.y) +
				(position.z - c->position.z) * (position.z - c->position.z);
			if (dist <= attractionRadius * c->getRadius() * attractionRadius * c->getRadius() && dist < minDist) {
				minDist = dist;
				closerCTL = c;
			}
		}
	}
	if (closerCTL != nullptr) {
		attractionDirection = (closerCTL->position - position);
		attractionDirection.normalize();
		state = attracted;
		return nullptr;
	} else {
		state = flying;
		return nullptr;
	}

	return nullptr;
}

VCell *CTLCell::attraction(double dt) {
	color[0] = 0.;
	color[1] = 0.7 * (1.0 - inhibitorySignal);
	color[2] = 0.7 * inhibitorySignal;
	uniform_real_distribution<double> dnoise(1.0, 1.0);
	setVelocity(attractionDirection * dt * attractionMoveSpeed * Scenario::MinutesPerTick);
	// looking for contact with tumor cells
	for (auto *c : getConnectedCells()) {
		if (c->type == Tumor) {
			state = contactWithTumor;
			tScan = 0;
			setVelocity(MecaCell::Vec(0,0,0));
			return nullptr;
		}
	}

	// looking for attractions
	attractionDirection.x = 0;
	attractionDirection.y = 0;
	attractionDirection.z = 0;
	double minDist = attractionRadius * getRadius() * attractionRadius * getRadius();
	CTLCell *closerCTL = nullptr;
	for (auto *c : scenario->collidingCTLCells) {
		if (c != this) {
			double dist = (position.x - c->position.x) * (position.x - c->position.x) + 
				(position.y - c->position.y) * (position.y - c->position.y) +
				(position.z - c->position.z) * (position.z - c->position.z);
			if (dist <= attractionRadius * c->getRadius() * attractionRadius * c->getRadius() && dist < minDist) {
				minDist = dist;
				closerCTL = c;
			}
		}
	}
	if (closerCTL != nullptr) {
		attractionDirection = (closerCTL->position - position);
		attractionDirection.normalize();
		state = attracted;
	} else {
		state = flying;
		return nullptr;
	}

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
				tScan += dt * Scenario::MinutesPerTick;
				// not yet detected => scanning
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

/*VCell *CTLCell::exhaust(double dt) {
	randomMove(dt);

	return nullptr;
}*/
	
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
	
	if (attractingCounter > 0) {
		attractingCounter -= attractingReduceSpeed * Scenario::MinutesPerTick * dt;
		if (attractingCounter <= 0.0 && (std::find(scenario->collidingCTLCells.begin(), scenario->collidingCTLCells.end(), this) != scenario->collidingCTLCells.end())) {
			scenario->collidingCTLCells.erase(std::remove(scenario->collidingCTLCells.begin(), scenario->collidingCTLCells.end(), this), scenario->collidingCTLCells.end());
		}
	}

	switch (state) {
		case flying:
			return fly(dt);
			break;
		case attracted:
			if (attractionRadius > 0.0001)return attraction(dt);
			else return fly(dt);
			break;
		case contactWithTumor:
			return touchingTumor(dt);
			break;
//		case exhausted:
//			return exhaust(dt);
//			break;
	}
	return nullptr;
}


