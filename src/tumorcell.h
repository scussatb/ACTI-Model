#ifndef __TUMORCELL_H__
#define __TUMORCELL_H__

#include <mecacell/mecacell.h>
#include <mecacell/basicworld.hpp>

#include <vcell.hpp>
#include "scenario.h"

class TumorCell : public VCell  {
	using VCell::VCell;
public:        
	Scenario* scenario;

	enum TumorState {
		cycling,
		contactWithCTL
	};

	enum TumorCyclePhase {
		G1S,
		G2M
	};

	enum TumorCellType {
		JY,
		D10
	};
	enum TumorSensitivityType {
		Sensitive,
		Persistant,
		Resistant
	};

	TumorState state;	
	static TumorCellType cellType;

	// Cycle variables
	double cycleDuration;
	double g1sDuration;
       	double g2mDuration;	
	static std::normal_distribution<double> cycleDurationNormDist;
	static double G1SRatio;
	double age;
	TumorCyclePhase cyclePhase;   

	// scanned state variables
	double scanDuration;
	static std::normal_distribution<double> scanDurationNormDist;

	// defending state variables
	double lifeLevel;
	double hitStrength;
	double repairLevel;
	static std::normal_distribution<double> repairLevelNormDist;
	static double inhibitoryStrength;
	static std::uniform_real_distribution<double> inhibitingDist;
	TumorSensitivityType sensitivity;

	TumorCell(MecaCell::Vec v);	

	~TumorCell();

	double getAdhesionWith(const VCell *c);
	
	VCell *updateBehavior(double dt); 

	VCell* isCycling(double dt);
       	VCell* isInContactWithCTL(double dt);
	void resetCycle();
	void init();
	void init(TumorCell* c);

	bool isDetected();
};

#endif
