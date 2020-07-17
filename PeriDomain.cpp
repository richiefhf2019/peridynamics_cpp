// Function Definitions

#include "PeriDomain.h"
#include "PeriMaterialModel.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <stdlib.h>

namespace periDynamics{

	//-------------------------------------------------------------------------
	// Default Constructor
	PeriDomain::PeriDomain() {
	    nPeriParticle = 0;
        TIMESTEP = 0.0;
    	printInterval = 0;
		periMaterial = new PeriMaterialModel();
	}
	
	// Overload Constructor
	PeriDomain::PeriDomain(int newnPeriParticle, REAL newTIMESTEP, int newprintInterval) {

		nPeriParticle = newnPeriParticle;
		TIMESTEP = newTIMESTEP;
		printInterval = newprintInterval;
		periMaterial = new PeriMaterialModel();
	}

	// Destructor
	PeriDomain::~PeriDomain() {
	
	 //   // free the spaces of these pointer vector
		//if(periParticleVec.size() > 0) {
		//	for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++) {
		//		delete (*pt);
		//	}
		//}
		//
		//if(bottomBoundaryVec.size() > 0) {
	 //   for(std::vector<PeriParticle*>::iterator pt=bottomBoundaryVec.begin(); pt!=bottomBoundaryVec.end(); pt++) {
		//delete (*pt);
	 //   }
		//}


		//if(cubicTopBoundaryVec.size() > 0) {
	 //   for(std::vector<PeriParticle*>::iterator pt=cubicTopBoundaryVec.begin(); pt!=cubicTopBoundaryVec.end(); pt++) {
		//delete (*pt);
	 //   }
		//}

	 //   delete periMaterial;
	
	 //   periParticleVec.clear();
	 //   bottomBoundaryVec.clear();
	 //   cubicTopBoundaryVec.clear();

	}

	//-------------------------------------------------------------------------
	// Accessor Functions
	int PeriDomain::getnPeriParticle() const {
		return nPeriParticle;
	}

	int PeriDomain::getPrintInterval() const {
		return printInterval;
	}

	REAL PeriDomain::getTIMESTEP() const {
		return TIMESTEP;
	}

	//-------------------------------------------------------------------------
	// Mutator Functions
	void PeriDomain::setnPeriParticle(int newnPeriParticle) {
		nPeriParticle = newnPeriParticle;
	}

	void PeriDomain::setPrintInterval(int newprintInterval) {
		printInterval = newprintInterval;
	}

	void PeriDomain::setTIMESTEP(REAL newTIMESTEP) {
		TIMESTEP = newTIMESTEP;
	}

	//-------------------------------------------------------------------------
	// Utility Functions

	void PeriDomain::initial(const char* inputFile){
		std::cout << "------------------------------------------------------------------------------" << std::endl;
		std::cout << "Problem Initilization " << std::endl;
		std::cout << "------------------------------------------------------------------------------" << std::endl;
		std::cout << "Read data file ..." << std::endl;
		readData(inputFile);
		std::cout << "Calculate particle volume ..." << std::endl;
		calcParticleVolume();
		// writeMeshCheckVolume("checkv.dat"); exit(1);
		std::cout << "Calculate horizon size ..." << std::endl;
		calcHorizonSize();
		std::cout << "Construct neighor list ..." << std::endl;
		constructNeighbor();
		std::cout << "Calculate Kinv ..." << std::endl;
		calcParticleKinv();
	    for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
			(*pt)->initial();
	    }

		// prescrible the essential boundary condition
		std::cout << "Prescribe the boundary condition ..." << std::endl;
		prescribeEssentialBoundaryCondition(0);

		// calculate the stress at each particle
		std::cout << "Calculate particle stress ..." << std::endl;
		calcParticleStress();

		// calculate the acceleration for each particle
		std::cout << "Calculate particle acceleration ..." << std::endl;
		calcParticleAcceleration();


		// traction boundary 
		ApplyExternalForce(0);

	    // apply initial velocity boundary condition, in this case give the cubicTopBoundaryVec particles initial velocities
	    //for(std::vector<PeriParticle*>::iterator pt=cubicTopBoundaryVec.begin(); pt!=cubicTopBoundaryVec.end(); pt++){
		//	(*pt)->setInitVelocity(Vec(0.0, 0.0, 1.0));
	    //}

	}

    void PeriDomain::prescribeEssentialBoundaryCondition(const int istep){

		for(std::vector<PeriParticle*>::iterator pt=bottomBoundaryVec.begin(); pt!=bottomBoundaryVec.end(); pt++){
			(*pt)->prescribeBottomDisplacement(0.0);	// fix z displacement in the z direction
		}

		//REAL dispz;
		//if(istep <= 200) {
		//	dispz = 1.8*0.05*double(istep)/200.;
		//}else {
		//	dispz = 1.8*0.05;
		//}

		//for(std::vector<PeriParticle*>::iterator pt=topBoundaryVec.begin(); pt!=topBoundaryVec.end(); pt++){
		//	(*pt)->prescribeTopDisplacement(dispz);	// fix z displacement in the z direction
		//}

    } // end perscribeEssentialBoundaryCondition()


    void PeriDomain::solve(const char* outputFile){
		// open the tecplot file for output
		std::ofstream ofs(outputFile);
		int iframe = 0;
		writeParticleTecplot(ofs,iframe);
		std::cout << "------------------------------------------------------------------------------" << std::endl;
		std::cout << "Start of the time loop " << std::endl;
		std::cout << "------------------------------------------------------------------------------" << std::endl;
		std::ofstream datafile("uxyz.dat");
		datafile.setf(std::ios::scientific, std::ios::floatfield);
		datafile.precision(10);
		datafile << "VARIABLES = \"Time step\", \"UX\", \"UY\", \"UZ\"" << std::endl;
		for(int istep = 1; istep <= nsteps; istep++) {

			
		    runFirstHalfStep();

		    prescribeEssentialBoundaryCondition(istep);

		    checkBondParticleAlive();

		    calcParticleStress();

		    calcParticleAcceleration();

		    ApplyExternalForce(istep);

		    runSecondHalfStep();

			if( istep % printInterval == 0) {
				std::cout << "*** current time step is	" << istep << std::endl;
				iframe++;
				writeParticleTecplot(ofs,iframe);
				datafile << istep
				<< std::setw(20) << periParticleVec[568]->getDisplacement().getx()
				<< std::setw(20) << periParticleVec[568]->getDisplacement().gety()
				<< std::setw(20) << periParticleVec[568]->getDisplacement().getz() << std::endl;
			}
			if( istep % 200 == 0) {
				writeDisplacementData("ux.dat","uy.dat","uz.dat");
			}

		} // time loop
		ofs.close();
		datafile.close();
		std::cout << "------------------------------------------------------------------------------" << std::endl;
		std::cout << "Simulation Finished !" << std::endl;
		std::cout << "------------------------------------------------------------------------------" << std::endl;
		writeDisplacementData("ux.dat","uy.dat","uz.dat");
    } // end solve()


	void PeriDomain::writeDisplacementData(const char *outputFilex, const char *outputFiley, const char *outputFilez) {
		// displacment along the x axis
		std::ofstream ofs(outputFilex);
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(10);
		ofs << "VARIABLES = \"X\", \"UX\"" << std::endl;
		for(int index = 0; index < 5; index++){
			int node = Uxindex[index];
			ofs << std::setw(20) << periParticleVec[node]->getInitPosition().getx()
				<< std::setw(20) << periParticleVec[node]->getDisplacement().getx() << std::endl;
		}
		ofs.flush();
		ofs.close();

		//dispalcement along the y axis
		ofs.open(outputFiley);
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(10);
		ofs << "VARIABLES = \"Y\", \"UY\"" << std::endl;
		for(int index = 0; index < 5; index++){
			int node = Uyindex[index];
			ofs << std::setw(20) << periParticleVec[node]->getInitPosition().gety()
				<< std::setw(20) << periParticleVec[node]->getDisplacement().gety() << std::endl;
		}
		ofs.flush();
		ofs.close();

		//dispalcement along the z axis
		ofs.open(outputFilez);
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(10);
		ofs << "VARIABLES = \"Z\", \"UZ\"" << std::endl;
		for(int index = 0; index < 39; index++){
			int node = Uzindex[index];
			ofs << std::setw(20) << periParticleVec[node]->getInitPosition().getz()
				<< std::setw(20) << periParticleVec[node]->getDisplacement().getz() << std::endl;
		}
		ofs.flush();
		ofs.close();


	}

    void PeriDomain::runFirstHalfStep(){

		for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
		    (*pt)->updateDisplacement(TIMESTEP);
		}

    } // end runFirstHalfStep()

    void PeriDomain::runSecondHalfStep(){

		for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
		    (*pt)->updateVelocity(TIMESTEP);
		}

    } // end runSecondHalfStep()


    void PeriDomain::constructNeighbor(){	
	// neighbor - searches and constructs the particle neighborlists
	// construct neighborlist for all particles ...
	// compute the weighting function for all particles ...

	for(std::vector<PeriParticle*>::iterator i_nt=periParticleVec.begin(); i_nt!=periParticleVec.end()-1; i_nt++){

	    Vec coord0_i = (*i_nt) ->getInitPosition();
	    REAL horizonSize_i = (*i_nt)->getHorizonSize();
	    for(std::vector<PeriParticle*>::iterator j_nt=i_nt+1; j_nt!=periParticleVec.end(); j_nt++){
			Vec coord0_j = (*j_nt)->getInitPosition();
			REAL tmp_length = vfabs(coord0_i-coord0_j);

			REAL horizonSize_j = (*j_nt)->getHorizonSize();
			REAL horizonSize_ij = (horizonSize_i+horizonSize_j)/2.0;//This will lead to the fact that horizion is not a sphere!!! 

			REAL ratio = tmp_length/horizonSize_ij;

			// establish the neighbor list
			if(ratio <= 2.0){ 

				// create bond
				Bond* bond_pt = new Bond(tmp_length, *i_nt, *j_nt);
				(*i_nt)->pushBackBondVec(bond_pt);
				(*j_nt)->pushBackBondVec(bond_pt);

				REAL factor = 3.0/(2.0*PI*horizonSize_ij*horizonSize_ij*horizonSize_ij); // for the factor of 3d window function
				
				// weighting function (influence function)
				if(ratio < 1.0){
					bond_pt->setWeight( factor*(2.0/3.0-ratio*ratio+0.5*ratio*ratio*ratio) );
				}
				else{
					bond_pt->setWeight( factor*(2.0-ratio)*(2.0-ratio)*(2.0-ratio)/6.0 );
				}		     
			} // if(ratio<2.0)

	    } // end j_nt
	} // end i_nt

    } // end constNeighbor()

    void PeriDomain::readData(const char *InputFile){
		// readData - reads controlling parameters, particle positions and mesh connectivities
		// @param char * - reference of the input file name
    	
		std::ifstream ifs(InputFile);
			if(!ifs) {
		    	std::cout << "stream error!" << std::endl; exit(-1);
			}
		
   		std::string tmp;
		// std::getline(ifs, tmp);	// read header in the first line
		
		ifs >> ndim >> nPeriParticle >> nele;
		
		// std::getline(ifs, tmp);	// read header in the third line
		// std::getline(ifs, tmp);	// read header in the forth line
		// std::getline(ifs, tmp);	// read header in the forth line
		
		ifs >> TIMESTEP >> nsteps >> printInterval >> rampStep;
		
		// std::getline(ifs, tmp);	// read header in the sixth line
		// std::getline(ifs, tmp);
		
		int typeConstitutive;
		ifs >> typeConstitutive;
		periMaterial->setTypeConstitutive(typeConstitutive);
		
		// std::getline(ifs, tmp);	// read header in the eighth line
		// std::getline(ifs, tmp);
		
		REAL periPoisson, periYoung, periDensity;
		ifs >> periPoisson >> periYoung >> periDensity >> bodyDensity;
		periMaterial->setPoisson(periPoisson);
		periMaterial->setYoung(periYoung);
		periMaterial->setDensity(periDensity);
		
		// std::getline(ifs, tmp);	// read header in the tenth line
		// std::getline(ifs, tmp);

		// read particle information, create and store PeriParticle objects into periParticleVec
		for(int ip=0; ip<nPeriParticle; ip++){
		    REAL tmp_x, tmp_y, tmp_z;
		    int tmp_int;
		    ifs >> tmp_int >> tmp_x >> tmp_y >> tmp_z;
		    PeriParticle* tmp_pt = new PeriParticle(tmp_x, tmp_y, tmp_z);
		    periParticleVec.push_back(tmp_pt);

		    // check bottom boundary particles
		    if( fabs(tmp_z+0.15)< 0.151 ){	// bottom particle
				bottomBoundaryVec.push_back(tmp_pt);
		    }
			else if( fabs(tmp_z-2.0)<0.01 ){	// top particle
				topBoundaryVec.push_back(tmp_pt);
		    }

		    //// check cubic top boundary particles
		    //if( fabs(tmp_z-0.013) < 1.0e-5 && tmp_x < 0.00375 && tmp_x > -0.00375 && tmp_y < 0.00375 && tmp_y > -0.00375){	// cubic top particle
			//cubicTopBoundaryVec.push_back(tmp_pt);
		    //}

		}
		
		// getline(ifs, tmp);	// read the header
		// getline(ifs, tmp);
		
		// read the connectivity information
		connectivity = new int*[nele];
		for(int iel=0; iel<nele; iel++){
			connectivity[iel] = new int[8];
		}
		
		for(int iel=0; iel<nele; iel++){
		    int tmp_int;
		    ifs >> tmp_int;
		    for(int node=0; node<8; node++){
		   	ifs >> connectivity[iel][node]; 
		    }
		}

		// read particle indices for linear elasticity verfication, can be deleted
		for(int ip=0; ip < 39; ip++) {ifs >> Uzindex[ip];}
		for(int ip=0; ip < 5; ip++) {ifs >> Uyindex[ip];}
		for(int ip=0; ip < 5; ip++) {ifs >> Uxindex[ip];}

		// calculate material parameters based on input
		periMaterial->calcMaterialParameters();

		setInitIsv();

    } // readData()

    void PeriDomain::writeMesh(const char *outputFile){
		std::ofstream ofs(outputFile);
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(10);
		ofs << "Title = \"Mesh Checking\"" << std::endl;
		ofs << "VARIABLES = \"X\", \"Y\",\"Z\"" << std::endl;
		ofs << "ZONE N = " << nPeriParticle << " E = " << nele << ", F = FEPOINT ET = BRICK" << std::endl;
		for(int node = 0; node < nPeriParticle; node++){
			ofs << std::setw(20) << periParticleVec[node]->getInitPosition().getx()
				<< std::setw(20) << periParticleVec[node]->getInitPosition().gety() 
				<< std::setw(20) << periParticleVec[node]->getInitPosition().getz() << std::endl;
		}
		for(int iel = 0; iel < nele; iel++){
			for(int node = 0; node < 8; node++){
				ofs << std::setw(10) << connectivity[iel][node]; 
			}
			ofs << std::endl;
		}
		ofs.close();
	} // end writeMesh


    void PeriDomain::writeMeshCheckVolume(const char *outputFile){
		std::ofstream ofs(outputFile);
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(10);
		ofs << "Title = \"Volume Checking\"" << std::endl;
		ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"V\"" << std::endl;
		ofs << "ZONE N = " << nPeriParticle << " E = " << nele << ", F = FEPOINT ET = BRICK" << std::endl;
		// Output the coordinates and the array information
		for(std::vector<PeriParticle*>::iterator pt = periParticleVec.begin(); pt!= periParticleVec.end(); pt++) {
			ofs << std::setw(20) << (*pt)->getInitPosition().getx()
				<< std::setw(20) << (*pt)->getInitPosition().gety() 
				<< std::setw(20) << (*pt)->getInitPosition().getz() 
				<< std::setw(20) << (*pt)->getParticleVolume()
				<< std::endl;
		}
		for(int iel = 0; iel < nele; iel++){
			for(int node = 0; node < 8; node++){
				ofs << std::setw(10) << connectivity[iel][node]; 
			}
			ofs << std::endl;
		}
		ofs.close();
	} // end writeMeshCheckVolume

    void PeriDomain::writeParticleTecplot(std::ofstream &ofs, const int iframe) {
		ofs.setf(std::ios::scientific, std::ios::floatfield);
		ofs.precision(10);
		if(iframe == 0) {
			ofs << "Title = \"Particle Information\"" << std::endl;
			ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" \"Vz\" \"KE\" \"P\" \"Mises\"" << std::endl;
		}
		ofs << "ZONE T =\" " << iframe << "-th Load Step\" "<< std::endl;
		// Output the coordinates and the array information
		REAL pressure, vonMisesStress;
		Matrix sigma;
		for(std::vector<PeriParticle*>::iterator pt = periParticleVec.begin(); pt!= periParticleVec.end(); pt++) {
			sigma = (*pt)->getSigma();
			pressure = sigma(1,1) + sigma(2,2) + sigma(3,3);
			vonMisesStress = sqrt(( (sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))
					         + (sigma(2,2)-sigma(3,3))*(sigma(2,2)-sigma(3,3))
							 + (sigma(1,1)-sigma(3,3))*(sigma(1,1)-sigma(3,3))
							 + sigma(1,2)*sigma(1,2)
							 + sigma(2,3)*sigma(2,3)
							 + sigma(3,1)*sigma(3,1) )/2.0);
			ofs << std::setw(20) << (*pt)->getInitPosition().getx() + (*pt)->getDisplacement().getx()
				<< std::setw(20) << (*pt)->getInitPosition().gety() + (*pt)->getDisplacement().gety()
				<< std::setw(20) << (*pt)->getInitPosition().getz() + (*pt)->getDisplacement().getz()
				<< std::setw(20) << (*pt)->getDisplacement().getx()
				<< std::setw(20) << (*pt)->getDisplacement().gety() 
				<< std::setw(20) << (*pt)->getDisplacement().getz() 
				<< std::setw(20) << (*pt)->getVelocity().getx()
				<< std::setw(20) << (*pt)->getVelocity().gety() 
				<< std::setw(20) << (*pt)->getVelocity().getz() 
				<< std::setw(20) << vfabs((*pt)->getVelocity())
				<< std::setw(20) << pressure
				<< std::setw(20) << vonMisesStress
				<< std::endl;
			ofs.flush();
		}
	}

    void PeriDomain::calcDeformationGradient(){	
		// calcDeformationGradient - calculates the deformation gradient for all peri-particles
	}

	void PeriDomain::calcHorizonSize(){

		for(int iel = 0; iel < nele; iel++){

			// get initial positions of the particles in this connectivity
 			int n1 = connectivity[iel][0]; 	// index of first node
			int n2 = connectivity[iel][1];
			int n4 = connectivity[iel][3];
			int n5 = connectivity[iel][4];

			Vec coord1 = periParticleVec[n1-1]->getInitPosition();	// the index of input file is starting from 1
			Vec coord2 = periParticleVec[n2-1]->getInitPosition();
			Vec coord4 = periParticleVec[n4-1]->getInitPosition();
			Vec coord5 = periParticleVec[n5-1]->getInitPosition();

			REAL tmp1, tmp2, tmp3;
			tmp1 = vfabs(coord2-coord1);
			tmp2 = vfabs(coord4-coord1);
			tmp3 = vfabs(coord5-coord1);

			REAL tmpmax;	// max number in tmp1, tmp2 and tmp3
			tmpmax = std::max( std::max(tmp1, tmp2), std::max(tmp2, tmp3) );

			for(int node = 0; node < 8; node++){
			    periParticleVec[connectivity[iel][node]-1]->replaceHorizonSizeIfLarger(1.5075*tmpmax);
			}

		}

	} // end calcHorizonSize()

	void PeriDomain::setInitIsv(){

	    REAL isv_tmp;
	    if(periMaterial->getTypeConstitutive() == 1){ // 1---implicit, 2---explicit
			isv_tmp = periMaterial->getChi();
	    }
	    else{
			isv_tmp = periMaterial->getC();
	    }

	    for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
			(*pt)->setInitIsv(isv_tmp);
	    }

       	} // setInitIsv()	

	void PeriDomain::calcParticleVolume() {
		
		int *numofpieces;
		REAL *particleVolume;
		REAL xi, eta, zeta;
		numofpieces = new int[nPeriParticle];
		for(int node = 0; node < nPeriParticle; numofpieces[node] = 0, node++);
		int nip = 2;
		Matrix gp_loc3D;
		Matrix gp_weight3D;
		gauss3D( nip, gp_loc3D, gp_weight3D);
		particleVolume = new REAL[nPeriParticle];
		for(int node = 0; node < nPeriParticle; particleVolume[node] = 0.0, node++);
		Matrix xl(3,8);
		Matrix shp;
		for(int iel = 0; iel < nele; iel++) {
			for(int node = 0; node < 8; node++) {
				int nposition = connectivity[iel][node]-1;
				xl(1,node+1) = periParticleVec[nposition]->getInitPosition().getx();
				xl(2,node+1) = periParticleVec[nposition]->getInitPosition().gety();
				xl(3,node+1) = periParticleVec[nposition]->getInitPosition().getz();
				numofpieces[nposition] += 1;
			}
			for(int ik = 0; ik < 8; ik++) {
				xi =   gp_loc3D(1,ik+1);
				eta =  gp_loc3D(2,ik+1);
				zeta = gp_loc3D(3,ik+1);
				// call fem function to get the volume
				REAL xsj = 0.0;
				shp3d(xi, eta, zeta, xl, shp, xsj);
				for(int node = 0; node < 8; node++) {
					int nposition = connectivity[iel][node] - 1;
					particleVolume[nposition] = particleVolume[nposition] + gp_weight3D(1,ik+1)*xsj/8.0;
				}
			}
		}

		//Commented to compare result with the fortran code
		for(int node = 0; node < nPeriParticle; node++) {
			particleVolume[node] = particleVolume[node]*8.0/(REAL(numofpieces[node]));
		}

		// store the particle volume into the object periParticle
		for(int node = 0; node < nPeriParticle; node++) {
			periParticleVec[node]->setParticleVolume(particleVolume[node]);
		}

		delete numofpieces;
		delete particleVolume;
	} // end calcParticleVolume

void PeriDomain::checkBondParticleAlive(){

	    // compute the bond length and check for whether a bond or particle is alive or not
	    for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){

		if( (*pt)->getIsAlive() ){ // particle not alive, then go to next particle
		    (*pt)->checkParticleAlive(periMaterial->getBondStretchLimit());
		}

	    } // end particle

} // end checkBondParticleAlive()

	
void PeriDomain::calcParticleKinv(){
	    // Compute the inverse of the shape tensor K
	    for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++) {
			(*pt)->calcParticleKinv();
	    } 

} // end calcParticleKinv()	



  void PeriDomain::calcParticleStress(){
	for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
	    (*pt)->calcParticleStress(periMaterial);
    	} 
   } // calcParticleStress()


    void PeriDomain::calcParticleAcceleration(){
		for(std::vector<PeriParticle*>::iterator pt=periParticleVec.begin(); pt!=periParticleVec.end(); pt++){
		    (*pt)->calcParticleAcceleration(periMaterial->getDensity());
		}

    } // calcParticleAcceleration()

	void PeriDomain::ApplyExternalForce(int istep) {
		// deal with the external force, applied at the top of the boundary
		REAL factor = 0.0;
		if(istep <= rampStep) {
			factor = REAL(istep)/REAL(rampStep);
		}
		else {
			factor = REAL(1.0);
		}

		for(std::vector<PeriParticle*>::iterator pt=topBoundaryVec.begin(); pt!=topBoundaryVec.end(); pt++){
			Vec newAccleration = (*pt)->getAcceleration() + Vec(0.0,0.0,factor*bodyDensity);
			(*pt)->setAcceleration(newAccleration);
		}
	}

	
} // end namespace periDynamics
