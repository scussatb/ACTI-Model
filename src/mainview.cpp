#include <mecacell/mecacell.h>
#include <cstring>
//#include <mecacell/viewer/viewer.h>
#include "scenario.h"
#include "ctlcell.h"
#include "tumorcell.h"
//#include <QtCore>
//#include <QtGui>

int main(int argc, char **argv) {
	Scenario::ctlTumorRatio = stod(argv[1]);
	CTLCell::hitStrength = 1.0 / stod(argv[2]);
	TumorCell::inhibitoryRadius = stod(argv[3]);
	CTLCell::inhibitoryProbability = stod(argv[4]);
	Scenario::simuDurationHour = stoi(argv[5]);
	TumorCell::cellType = std::strcmp(argv[6],"JY")==0?TumorCell::JY:TumorCell::D10;
	TumorCell::cycleDurationNormDist=std::normal_distribution<double>(TumorCell::cellType==TumorCell::D10?1200:1065, 60);

	if (argc==7) {
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
