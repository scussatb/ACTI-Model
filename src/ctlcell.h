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
		contactWithTumor,
	};

	CTLState state;

	// Attacking variables
	double inhibitorySignal;
	static double inhibitoryProbability;
	static double hitStrength;

	// Flying variable
	static double flySpeed;
	static std::exponential_distribution<double> flyExpDist;
	static std::normal_distribution<double> flyBrownianDistanceNormDist;
	static std::uniform_real_distribution<double> flyAngleUniDist;
	static pareto_distribution flyLengthParDist;
	MecaCell::Vec flyDirection;
	double flyDuration;
	double nextChange;
	int nChange;
	
	int nbKilledCell;
	int tumorSynapticTime;
	int ctlSynapticTime;

	void randomMove(double dt);

	VCell *fly(double dt);
	VCell *touchingTumor(double dt);

	CTLCell(MecaCell::Vec v);	

	~CTLCell();

	double getAdhesionWith(const VCell *c);
	
	VCell *updateBehavior(double dt); 
};

#endif
