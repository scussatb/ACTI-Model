#ifndef __SCENARIO_H__
#define __SCENARIO_H__

#include <mecacell/mecacell.h>
//#include <mecacell/viewer/viewer.h>
#include "vcell.hpp"
//#include "tumorcell.h"
//#include "ctlcell.h"

#define BOXSIZE 375

using World = MecaCell::BasicWorld<VCell, MecaCell::Euler>;

class TumorCell;
class CTLCell;

class Scenario {

private:
	World w;

public:
//	static MecacellViewer::Viewer<Scenario> viewer;
	static double BOX_COEF;
	static double BOX_SIZE;
	static double TIME_ACCEL;
	static double MinutesPerTick;
	static double HoursPerTick;

	// stat options
	static bool synapticContacts;
	static bool hourlyReports;
	static bool finalReport;
	static bool perCapitaKilling;
	
	static int reinjectAfter;

	int currentHour;
	int oldHour;
	int delay;

	int nTumorCells;
	static double ctlTumorRatio;
	std::vector<TumorCell*> tumorCells;
	unsigned int nDeadTumorCells;
	
	static int simuDurationHour;
	
	std::vector<CTLCell*> ctlCells;
	std::vector<CTLCell*> collidingCTLCells;
	int nInhibitedCells;

	Scenario(int nTumorCells);
	Scenario() {nTumorCells = 1000;};
	~Scenario() {
		tumorCells.clear();
		collidingCTLCells.clear();
		ctlCells.clear();
	};

	void init(int argc, char **argv);
	void loop();
	World &getWorld();
};

#endif
