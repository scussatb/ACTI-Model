#include "tumorcell.h"
#include "ctlcell.h"

#include <limits>


TumorCell::TumorCellType TumorCell::cellType = D10;
std::normal_distribution<double> TumorCell::cycleDurationNormDist(TumorCell::cellType==D10?1200:1065, 60); // cycles in 20h/17.75h +- 1h
double TumorCell::G1SRatio = TumorCell::cellType==D10?0.70:0.77;
double TumorCell::inhibitoryRadius = 87.5;
std::uniform_real_distribution<double> TumorCell::inhibitingDistribution(0, 1);


TumorCell::TumorCell(MecaCell::Vec v):	
	VCell(v) {
	type = Tumor;
}

void TumorCell::init() {
	lifeLevel = 1.0;
	resetCycle();
	type = Tumor;
	color[0] = 0.0;
	color[1] = 0.0;
	color[2] = 0.0;
}


void TumorCell::init(TumorCell *c) {
	scenario = c->scenario;
	w = c->w;
	lifeLevel = c->lifeLevel;
	resetCycle();
	type = Tumor;
}	

TumorCell::~TumorCell() {
}

double TumorCell::getAdhesionWith(const VCell *c) {
	return (c->type == Tumor)?0.8:0.8;
}	
	
VCell *TumorCell::updateBehavior(double dt) {
	if (isDead()) {
		return nullptr;
	}
	if (position.x * position.x + position.z * position.z > 0.4*0.4*Scenario::BOX_SIZE * Scenario::BOX_SIZE/4.) {
		scenario->tumorCells.erase(std::remove(scenario->tumorCells.begin(), scenario->tumorCells.end(), this), scenario->tumorCells.end());
		die();
		return nullptr;
	}
	if (getVelocity().y<-10.0) {
		// currently falling, do nothing
		return nullptr;
	}
	
	switch (state) {
		case cycling:
			return isCycling(dt);
			break;
		case contactWithCTL:
			return isInContactWithCTL(dt);
			break;
	}
	return nullptr;
}

void TumorCell::resetCycle() {
	age = 0.;
	cycleDuration = cycleDurationNormDist(MecaCell::globalRand);
	g1sDuration = cycleDuration * G1SRatio;
	g2mDuration = cycleDuration - g1sDuration;
	cyclePhase = G1S;
	state = cycling;
}

VCell* TumorCell::isCycling(double dt) {
	color[0] = 0.0;
	color[1] = 0.0;
	color[2] = 0.0;
	
	age += Scenario::MinutesPerTick * dt;
	if (age > g1sDuration) {
		// currently in G2M
	 	if (getRelativeVolume() < 2.0) {
			color[0] = getRelativeVolume()-1.0;
			color[1] = getRelativeVolume()-1.0;
			color[2] = getRelativeVolume()-1.0;
			// cell is in G2M and no big enough => growing
			grow(dt / g2mDuration * Scenario::MinutesPerTick);
	       } else {
			//cell is big enough => divide
			TumorCell* c = divide<TumorCell>();
			c->init(this);
			resetCycle();
			scenario->tumorCells.push_back(c);
			return c;
	       }
	}

	// contact with CTL
	for (auto *c : getConnectedCells()) {
		if (c->type == CTL) {
			state = contactWithCTL;
			return nullptr;
		}
	}

	return nullptr;
}

VCell* TumorCell::isInContactWithCTL(double dt) {
	color[0] = 0.0;
	color[1] = 0.0;
	color[2] = 0.0;

	int nCTLConnected = 0;
	for (auto *c : getConnectedCells()) {
		if (c->type == CTL) {
			CTLCell *ctl = (CTLCell*) c;
			nCTLConnected++;
			// is attacked
			if (lifeLevel > 0) {
				if (ctl->inhibitorySignal<=0.01) {
					lifeLevel -= dt * Scenario::MinutesPerTick * ctl->hitStrength;
					color[0] = 0.7;
					color[1] = 0.0;
					color[2] = 0.0;
					for (auto* cs : scenario->ctlCells) {
						if (cs != c) {
							if (cs->inhibitorySignal < 0.01) {
								double dist = 	std::sqrt((cs->getPosition().x - position.x) * (cs->getPosition().x - position.x) +
							     		(cs->getPosition().y - position.y) * (cs->getPosition().y - position.y) +
								(cs->getPosition().z - position.z) * (cs->getPosition().z - position.z));
								if (dist<inhibitoryRadius) {
									if (inhibitingDistribution(MecaCell::globalRand) < CTLCell::inhibitoryProbability) {
										cs->inhibitorySignal = 1.0; 
										scenario->nInhibitedCells++;
									}	
								}
							}
						}
					}
				}
			} else {
				scenario->nDeadTumorCells++;
				scenario->tumorCells.erase(std::remove(scenario->tumorCells.begin(), scenario->tumorCells.end(), this), scenario->tumorCells.end());
				die();
				ctl->nbKilledCell++;
				break;
			}
		}
	}

	// test if still attacked/scanned
	if (nCTLConnected <= 0) {
		state = cycling;
		return nullptr;
	}

	return nullptr;
}

bool TumorCell::isDetected() {
	return true;
}
