#ifndef __VCELL_HPP__
#define __VCELL_HPP__

#include <mecacell/mecacell.h>
#include <mecacell/basicworld.hpp>

class TumorCell;
class CTLCell;

class VCell : public MecaCell::ConnectableCell<VCell> {
	using Base = MecaCell::ConnectableCell<VCell>;
	using Base::Base;
public:        
	enum CellType {
		Tumor,
		CTL
	};


	VCell (MecaCell::Vec v) :
		ConnectableCell::ConnectableCell(v) {};

	// shared variables
	double age;
	MecaCell::BasicWorld<VCell, MecaCell::Euler> *w;
	CellType type;

	virtual ~VCell() {};

	virtual double getAdhesionWith(const VCell *c) = 0;
	
	virtual VCell *updateBehavior(double dt) = 0; 
};



#endif
