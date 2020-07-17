// Header ==> Function Declarations
#pragma once
#include <iostream>
#include <string>
#include <vector>

#include "globfuncs.h"
#include "realtypes.h"
#include "PeriParticle.h"
#include "PeriMaterialModel.h"
#include "Matrix.h"


#ifndef PERIDOMAIN_H
#define PERIDOMAIN_H

namespace periDynamics{

// PeriDomain is the assembly of particles, it's very similar to the class assembly in dem
// all the I(nitialize)R(un)F(inalize) are done in this class
class PeriDomain {

public:
	//-------------------------------------------------------------------------
	//Default Constructor
    PeriDomain();

	//Overload Constructor
	PeriDomain(int, REAL, int);

	//Destructor
	~PeriDomain();

	//-------------------------------------------------------------------------
	//Accessor Functions
	int getnPeriParticle() const;
		// getnPeriParticle - returns number of peri-particles
		// @return - number of peri-particles

	REAL getTIMESTEP() const;
		// getTIMESTEP - returns the simulation timestep
		// @return - the simulation timestep

	int getPrintInterval() const;
		// getprintInterval - returns the output frequency
		// @return - the output frequency
	
	//-------------------------------------------------------------------------
	// Mutator Functions
	void setnPeriParticle(int);
		// setnPeriParticle - sets the number of peri-particles
		// @param int - number of peri-particles

	void setTIMESTEP(REAL);
		// setTIMESTEP - sets the time step of the simulation
		// @param REAL - time step of the simulation

	void setPrintInterval(int);
		// setprintInterval - sets the output frequency of the simulation
		// @param int - output frequency of the simulation

	void setInitIsv();		
		
	//-------------------------------------------------------------------------
	// Utility Functions
	void initial(const char* );	
		// initial - initializes the velocity, displacement and acceleration 
		// @param const char* - name of the input file 

    void prescribeEssentialBoundaryCondition(const int);
		// apply the fixed boundary condition on the bottom particles

    void solve(const char*);	
		// solve - solve this problem
		// @param const char* - output file name for tecplot visualization

    void runFirstHalfStep();
		// run step1 and step2 in Page 5 of Houfu's notes, equations (15) and (16)

    void runSecondHalfStep();
		// run step3 in page 5 of Houfu's notes, equation (17)

    void constructNeighbor();	
		// neibhor - searches and constructs the particle neighborlists

    void readData(const char *);	
		// readData - reads controlling parameters, particle positions and mesh connectivities
		// @param char * - reference of the input file name

    void writeMesh(const char *);
		// writeMesh - outputs the mesh, used for problem checking
		// @param char * - reference of the output file name 

    void writeMeshCheckVolume(const char *);
		// writeMeshCheckVolume - outputs the mesh and volume, will be used for volume comuptation checking
		// @param char * - reference of the output file name 

	void writeParticleTecplot(std::ofstream &, const int);
		// writeParticleTecplot - outputs the information of all particles, for tecplot visualization
		// @param (std::ofstream, int) - (reference of the output file name, frame index)

    void calcDeformationGradient();	
		// calcDeformationGradient - calculates the deformation gradient for each peri-particle

	void calcHorizonSize();	
		// calcHorizonSize - calculates the horizon size for each peri-particle

	void calcParticleVolume();
		// calcParticleVolume - calculates the particle volume associated with each particle

    void calcParticleKinv();	
		// calculate the inverse of K for all particles	

    void calcParticleStress();
		// calculate the Cauchy Stress for all particles

    void calcParticleAcceleration();
		// calculate the acceleration for all particles		

    void checkBondParticleAlive();
		// for each particle, check the state(alive or not) of each surrounding bond;
		// if there's no alive bond for a particle, then this particle is disabled.
	void ApplyExternalForce(int istep);

	// for elasticity verification purpose.
	int Uzindex[39]; // particle index along the z direction.[x = 0, y = 0]
	int Uyindex[5];  // particle index along the y direction.[x = 0, z = 0]
	int Uxindex[5];  // particle index along the x direction.[y = 0, z = 0]

	void writeDisplacementData(const char *, const char *, const char *);
 
private:
	// Member Variables
    int nPeriParticle;	// number of PeriParticles in the domain
	int nele; // number of elements in the mesh

    REAL TIMESTEP;	// simulation time step

    int ndim;	// dimension of the problem
    int nsteps;	// number of total steps
	int rampStep; // steps for the applied boundary to reach maximum
    int printInterval;	// print interval
    int **connectivity;	// mesh connectivity

	REAL bodyDensity; // force density applied at the boundary

    std::vector<PeriParticle*> periParticleVec;	// coordinates of all the particles in this domain
    std::vector<PeriParticle*> bottomBoundaryVec;	// particles that are in the bottom boundary
	std::vector<PeriParticle*> topBoundaryVec;	// particles that are in the bottom boundary
    std::vector<PeriParticle*> cubicTopBoundaryVec;	// particles that are in the top boundary of the cubic
//    std::vector<Bond*> totalBondVec;	// all the Bonds in the domain

    PeriMaterialModel* periMaterial;	

}; // end PeriDomain


} // end periDynamics

#endif
