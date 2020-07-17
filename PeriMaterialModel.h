#pragma once
#ifndef PERIMATERIALMODEL_H
#define PERIMATERIALMODEL_H

#include "realtypes.h"
#include "parameter.h"
#include "Matrix.h"


namespace periDynamics{


class PeriMaterialModel {

public:
    PeriMaterialModel() {
	Poisson = 0;
	Young = 0;
	lambda = 0;
	kBulk = 0;
	tangentModulus = zeros(6,6);

	hchi = 0;
	chi = 0;
 	c = 0;
	phi = 0;
	psi = 0;
	kappa = 0;
	rEllip = 0;
	beta = 0;
	bondStretchLimit = 0;
	typeConstitutive = 1;
	density = 0;

    }

    REAL getPoisson() const {return Poisson;}
    REAL getYoung() const {return Young;}
    REAL getDensity() const {return density;}
    REAL getLambda() const {return lambda;}
    REAL getChi() const {return chi;}
    REAL getC() const {return c;}
	REAL getBondStretchLimit() const {return bondStretchLimit;}
    REAL getAphi() const {return Aphi;}
    REAL getBphi() const {return Bphi;}
    REAL getApsi() const {return Apsi;}
    REAL getBpsi() const {return Bpsi;}
    REAL getMu() const {return mu;}
    REAL getKBulk() const {return kBulk;}
    REAL getHchi() const {return hchi;}
    int getTypeConstitutive() const {return typeConstitutive;}
    Matrix getTangentModulus() const {return tangentModulus;}

    void setPoisson(REAL poi) {Poisson = poi;}
    void setYoung(REAL you) {Young = you;}
    void setDensity(REAL den) {density = den;}
    void setTypeConstitutive(int type) {typeConstitutive = type;}
    void setBondStretchLimit(REAL limit) {bondStretchLimit = limit;}

    void calcMaterialParameters();



private:

    REAL Poisson;	// Poisson's ratio in peri-domain
    REAL Young;		// Young's modulus in peri-domain

    REAL lambda;
    REAL mu;
    REAL kBulk;
    
    Matrix tangentModulus;	// tangent modulus
    REAL hchi;
    
    REAL chi;	// kPa
    REAL c;	// kPa
    REAL phi;	// radians
    REAL psi;  // radians

    REAL kappa;// kPa
    REAL rEllip;// dimensionless

    REAL beta;	// TC

    REAL bondStretchLimit;

    REAL Aphi;
    REAL Bphi;
    REAL Apsi;
    REAL Bpsi;

    REAL density;

    int typeConstitutive;	// the type of constitutive relationship

};

} // end periDynamics

#endif
