#include "scenario.h"
#include <chrono>
#include <random>

#include "tumorcell.h"
#include "ctlcell.h"

double Scenario::BOX_COEF = 300.;
double Scenario::BOX_SIZE = 6.35 * BOX_COEF;
double Scenario::TIME_ACCEL = 1000;
double Scenario::MinutesPerTick = 60;
double Scenario::ctlTumorRatio = 1.0;
int Scenario::simuDurationHour = 18;

bool Scenario::synapticContacts = false;
bool Scenario::hourlyReports = true;
bool Scenario::finalReport = false;
bool Scenario::perCapitaKilling = false;

//MecacellViewer::Viewer<Scenario> Scenario::viewer;

Scenario::Scenario(int nTumorCells) {
	this->nTumorCells = nTumorCells;
}

void Scenario::init(int argc, char **argv) {
	currentHour = 0;
	oldHour = 0;
	delay = 0;
	tumorCells.clear();
	ctlCells.clear();
	nDeadTumorCells = 0;
	collidingCTLCells.clear();
	nInhibitedCells = 0;

	w.setG(MecaCell::Vec(0.0, -5.1, 0.0));
	w.setViscosityCoef(0.001);

	string plateFileName("../UBottomPlate_blender.obj");
	w.addModel("petri_dish", plateFileName);
	w.models.at("petri_dish").scale(MecaCell::Vec(BOX_COEF, 0.33*BOX_COEF, BOX_COEF));
	w.models.at("petri_dish").translate(MecaCell::Vec(0, -3*BOX_COEF, 0));

	std::uniform_real_distribution<double> d(-BOX_SIZE/2., BOX_SIZE/2.);
	int nCTLs = nTumorCells * ctlTumorRatio;
	double tumorCellRadius = 7.5;
	double ctlCellRadius = 3.0;
	double squaredInternalRadius = tumorCellRadius * tumorCellRadius * nTumorCells + ctlCellRadius * ctlCellRadius * nCTLs; 
	squaredInternalRadius *= 0.75;

	for (int i=0; i<nTumorCells; i++) {
		double x = d(MecaCell::globalRand);
		double z = d(MecaCell::globalRand);
		double y = BOX_COEF + d(MecaCell::globalRand) / 10.0; 
		while (x*x+z*z > squaredInternalRadius) {
			// outside the plate
			x = d(MecaCell::globalRand);
			z = d(MecaCell::globalRand);
		}
		TumorCell *c = new TumorCell(MecaCell::Vec(x, y, z));
        c->w = &w;
		c->scenario = this;
		c->setBaseRadius(tumorCellRadius);
		c->setRadius(tumorCellRadius);
		c->setMass(10);
		c->type = VCell::Tumor;
		c->init();
		c->setVelocity(MecaCell::Vec(0.0, -100.0, 0.0));
		std::exponential_distribution<double> ageDist(std::log(2.0));
		double tirage = ageDist(MecaCell::globalRand);
		while (tirage > 1.0) tirage = ageDist(MecaCell::globalRand);
		c->age = tirage * c->cycleDuration;
		if (c->age > c->g1sDuration) {
			c->grow((c->age-c->g1sDuration)/c->g2mDuration);
		}
		w.addCell(c);
		tumorCells.push_back(c);
	}
	squaredInternalRadius *= 1.1;
	for (int i=0; i<nCTLs; i++) {
		double x = d(MecaCell::globalRand);
		double z = d(MecaCell::globalRand);
		while (x*x+z*z > squaredInternalRadius) {
			// outside the plate
			x = d(MecaCell::globalRand);
			z = d(MecaCell::globalRand);
		}
		double y = BOX_COEF + d(MecaCell::globalRand) / 10.0; 
		CTLCell *c = new CTLCell(MecaCell::Vec(x, y, z));
	    c->w = &w;
		c->scenario = this;
		c->setBaseRadius(ctlCellRadius);
		c->setRadius(ctlCellRadius);
		c->setVelocity(MecaCell::Vec(0.0, -100.0, 0.0));
		w.addCell(c);
		ctlCells.push_back(c);
	}

}

void Scenario::loop() {
	// update tumor center of mass
	w.update();
	if (nDeadTumorCells > 0) {
		currentHour = (w.getNbUpdates()-delay) * w.getDt() * MinutesPerTick / 60;
	} else {
		delay++;
	}
		
	if (currentHour != oldHour) {
		oldHour = currentHour;
		if (hourlyReports) {
			cout << currentHour << "\t" << tumorCells.size() << "\t" << nDeadTumorCells << "\t" << (double)nDeadTumorCells / (tumorCells.size()+nDeadTumorCells) << "\t" << ctlCells.size() << "\t" << (double)nDeadTumorCells/ctlCells.size() << "\t" << nInhibitedCells << "\t" << (double)nInhibitedCells/ctlCells.size() << endl;
		}
	}

}

World &Scenario::getWorld() { return w; }
