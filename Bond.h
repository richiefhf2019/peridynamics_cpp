// Header ==> Function Declarations
#include<iostream>
#include<string>

//#include "PeriParticle.h"
#include "realtypes.h"
#include "Vec.h"
#include "Matrix.h"

#ifndef BOND_H
#define BOND_H

namespace periDynamics {

class PeriParticle;

class Bond {

public:
	//-------------------------------------------------------------------------
	// Default Constructor
	Bond();
	
	// Overload Constructor
	Bond(REAL , PeriParticle* , PeriParticle* );
	
	// Destructor
	~Bond();

	//-------------------------------------------------------------------------
	// Accessor Functions
	bool getIsAlive() const;
		// getisAlive - returns state of the Bond
		// @return - state of the Bond

	REAL getWeight() const;
		// getweight - returns weight of the Bond
		// @return - weight of the Bond

	REAL getInitLength() const;
		// getinitLength - returns initial length of the Bond
		// @return - initial length of the Bond
		
	REAL getParticleVolume(bool) const;

	Vec getXi(bool) const;
	
 	Vec getEta(bool) const;

	Vec getEtaHalf(bool, const REAL) const;	

	PeriParticle* getPt1() const {return pt1;}

	PeriParticle* getPt2() const {return pt2;}
 
	Matrix getMicroK(const bool) const;	// get the contribution of K from one single bond

	Matrix getMicroN(const bool, const bool) const;

	Matrix getMicroNHalf(const bool, const bool, const REAL) const;

	Matrix getMicroNDeltaU(const bool, const bool, const REAL) const;

	//-------------------------------------------------------------------------
	// Mutator Functions
	void setIsAlive(bool) ;
		// setisAlive - sets state of the Bond
		// @param - state of the Bond

	void setWeight(REAL) ;
		// setweight - sets weight of the Bond
		// @param - weight of the Bond

	void setInitLength(REAL);
		// setinitLength - sets initial length of the Bond
		// @param - initial length of the Bond

	void setAliveFalse() {isAlive = false;}		
		
	//-------------------------------------------------------------------------
	// Utility Functions

	REAL calcCurrentLength();

private:
	// Member Variables
	bool isAlive; // if the Bond is alive or not
	REAL weight; // influence function
	REAL initLength; // initial Bond length
	PeriParticle *pt1; // what for? store address of the particles it belongs to.
	PeriParticle *pt2; 
};  // end Bond

}   // end periDynamics
#endif

