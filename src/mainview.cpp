#include <mecacell/mecacell.h>
#include <cstring>
//#include <mecacell/viewer/viewer.h>
#include "scenario.h"
#include "ctlcell.h"
#include "tumorcell.h"
//#include <QtCore>
//#include <QtGui>

int main(int argc, char **argv) {
	CTLCell::attractionRadius = stod(argv[1]);
	CTLCell::attractingReduceSpeed = 1.0 / stod(argv[2]);
	Scenario::ctlTumorRatio = stod(argv[3]);
	CTLCell::hitStrengthNormDist = std::normal_distribution<double>(1.0/stod(argv[4]), stod(argv[5]));
	TumorCell::repairLevelNormDist = std::normal_distribution<double>(stod(argv[6]), stod(argv[7]));
	TumorCell::inhibitoryStrength = stod(argv[8]);
	CTLCell::inhibitoryRecoverySpeed = stod(argv[9]);
	Scenario::simuDurationHour = stoi(argv[10]);
	TumorCell::cellType = std::strcmp(argv[11],"JY")==0?TumorCell::JY:TumorCell::D10;
	TumorCell::cycleDurationNormDist=std::normal_distribution<double>(TumorCell::cellType==TumorCell::D10?1200:1065, 60);
	Scenario::reinjectAfter = stoi(argv[12]);

	if (argc==13) {
		// console mode
		Scenario s(1000);
		s.init(argc, argv);
		double nInhibitedCTL = 0;
		while (s.currentHour != Scenario::simuDurationHour) {
			s.loop();
		}
		exit(0);
	}/* else {
		// with GUI
		MecacellViewer::Viewer<Scenario> v;
		return Scenario::viewer.exec(argc, argv);
	}*/
}
