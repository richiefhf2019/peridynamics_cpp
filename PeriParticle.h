//  Header ==> Function Declarations
#pragma once
#include <iostream>
#include <vector>

#include "Bond.h"
#include "Vec.h"
#include "Matrix.h"
#include "realtypes.h"
#include "globfuncs.h"
#include "PeriMaterialModel.h"

#ifndef PERIPARTICLE_H
#define PERIPARTICLE_H

namespace periDynamics{

class PeriDomain;

class PeriParticle {

public:
	// Default Constructor
    PeriParticle();
    PeriParticle(REAL x, REAL y, REAL z);

	~PeriParticle();
    
    //void calcKinv();	// calculate Matrix the inverse of K Matrix
    //void calcAcceleration();	// calculate the acceleration of the particle
    //void calcStress();		// calculate stress 

	void setParticleVolume(REAL newParticleVolume);
		// setParticleVolume - sets the volume of the particle
		// @param newParticleVolume - volume of the particle
		
    REAL getParticleVolume() const {return particleVolume;}
	REAL getHorizonSize() const {return horizonSize;}
    Vec getInitPosition() const {return initPosition;}
	Vec getDisplacement() const {return displacement;}
	Vec getVelocity() const {return velocity;}
	bool getIsAlive() const {return isAlive;}
	Matrix getSigma() const {return sigma;}
    Matrix getDeformationGradient() const {return deformationGradient;}	
    Matrix getParticleKinv() const {return Kinv;}
    Matrix getIsv() const {return isv;}
	Vec getVelocityHalf() const {return velocityHalf;}
	
	Vec getAcceleration() const {return acceleration;}
	
	
    void checkParticleAlive(REAL);
    void replaceHorizonSizeIfLarger(REAL tmp);	// replace this->horizonSize if tmp is larger
    void prescribeBottomDisplacement(REAL disp ) {displacement.setz(disp); velocity.setz(0.0);}
	void prescribeTopDisplacement(REAL disp) {displacement.setz(disp);}
    void setInitVelocity(const Vec& tmp) {velocity = tmp;}
    void setInitIsv(REAL isv_tmp) {isv(1,1) = isv_tmp;}

	void setAcceleration(Vec newAcceleration) {acceleration = newAcceleration;}
//    void pushBackNeighborVec(PeriParticle* pt) {neighborVec.push_back(pt);}
    void pushBackBondVec(Bond* bt) {bondVec.push_back(bt);}
	void setAliveFalse() {isAlive = false;}
	void setStressZero() {sigma = zeros(3,3);}
	void setAccelerationZero() {acceleration = 0;}
    void setParticleKinv(Matrix& tmp) {Kinv = tmp;}
    void calcParticleKinv();
    void calcParticleStress(const PeriMaterialModel*);
    void calcParticleAcceleration(const REAL);

    void updateDisplacement(const REAL);	// update displacement and velocityHalf
    void updateVelocity(const REAL);		// update velocity
    void initial();	
	// initial displacement, velocity and acceleration

	
private:
    bool isAlive; 	// if the peri-particle is alive

    //REAL particleDensity;	// particle density ==> defined globally, same material, usually
    //dem::Vec currPosition;	// current position vector of particle, no need, save spaces
    Vec initPosition;	// initial position vector of particle

    REAL particleVolume; // particle volume
    Vec displacement; // particle displacement
    Vec velocity;     // particle velocity
	Vec velocityHalf; // velocity at half time-step, for Velocity-Verlet integration
    Vec acceleration; // particle acceleration

    Matrix sigma;		// Cauchy stress
    Matrix deformationGradient;// deformation gradient tensor
    Matrix deformationGradientHalf;
    Matrix Kinv;	// the inverse K Matrix
    Matrix isv;		// a 1x5 matrix 

	REAL horizonSize;	// used to get neighbor list

//    std::vector<PeriParticle*> neighborVec;	// neighbor list of this particle
    std::vector<Bond*> bondVec;	// Bonds connected to this particle


}; // end particle


} // end periDynamics

#endif
