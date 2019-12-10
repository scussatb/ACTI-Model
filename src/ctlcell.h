#ifndef __CTLCELL_H__
#define __CTLCELL_H__

#include <mecacell/mecacell.h>
#include <mecacell/basicworld.hpp>
#include "paretodist.hpp"
#include "scenario.h"

#include "vcell.hpp"

class CTLCell : public VCell {
public:        
	Scenario* scenario;

	enum CTLState {
		flying,
		attracted,
		contactWithTumor,
		exhausted
	};

	enum CTLFlyType {
		brownian,
		jump
	};

	CTLState state;

	// Attracted state variables
	static double attractionRadius;
	MecaCell::Vec attractionDirection;
	static double attractionMoveSpeed;
	static double attractingReduceSpeed;
	double attractingCounter;
	// Attacking variables
	double tScan;
	static double tScanMax;
	double inhibitorySignal;
	static double inhibitoryRecoverySpeed;
	double hitStrength;
	static std::normal_distribution<double> hitStrengthNormDist;
	// Flying variable
	double flyTimeMode;
	static std::exponential_distribution<double> flyExpDist;
	static std::normal_distribution<double> flyBrownianDistanceNormDist;
	static std::uniform_real_distribution<double> flyAngleUniDist;
	static pareto_distribution flyLengthParDist;
	static std::exponential_distribution<double> flyJumpDurationExpDist;
	CTLFlyType flyType;
	MecaCell::Vec flyDirection;
	double flyDuration;
	double nextChange;
	int nChange;
	
	int nbKilledCell;
	int tumorSynapticTime;
	int ctlSynapticTime;

	void randomMove(double dt);

	VCell *fly(double dt);
	VCell *attraction(double dt);
	VCell *touchingTumor(double dt);
//	VCell *exhaust(double dt);

	CTLCell(MecaCell::Vec v);	

	~CTLCell();

	double getAdhesionWith(const VCell *c);
	
	VCell *updateBehavior(double dt); 
};

#endif
