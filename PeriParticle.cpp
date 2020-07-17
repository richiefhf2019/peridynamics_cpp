// Function Definitions
#include "PeriParticle.h"

namespace periDynamics {

	PeriParticle::PeriParticle(REAL x, REAL y, REAL z){
	    isAlive = true; 
	    initPosition.setx(x); initPosition.sety(y); initPosition.setz(z);
	    particleVolume = 0.0; 
	    displacement = 0.0;
	    velocity = 0.0;     
	    velocityHalf = 0.0; 
	    acceleration = 0.0; 
	
	    sigma = zeros(3,3);
	    deformationGradient = zeros(3,3);
	    deformationGradientHalf = zeros(3,3);
	    Kinv = zeros(3,3);
	    isv = zeros(1,5);
	
	} // end PeriParticle()
	
	
	PeriParticle::~PeriParticle(){
	
	    // free the spaces of these pointer vector
//	    for(std::vector<PeriParticle*>::iterator ip=neighborVec.begin(); ip!=neighborVec.end(); ip++){
//		delete (*ip);
//	    }
	    for(std::vector<Bond*>::iterator ib=bondVec.begin(); ib!=bondVec.end(); ib++) {
		delete (*ib);
	    }
	
//	    neighborVec.clear();
	    bondVec.clear();	
	
	} // end PeriParticle()
	
	
	void PeriParticle::setParticleVolume(REAL newParticleVolume) {
	
		particleVolume = newParticleVolume;
	
	} // end setParticleVolume
	
	
	void PeriParticle::replaceHorizonSizeIfLarger(REAL tmp) {
	
	    if(horizonSize < tmp){
			horizonSize = tmp;
	    }
	
	} // end replaceHorizonSizeIfLarger()


	void PeriParticle::calcParticleKinv(){

	    Matrix K(3,3);
		for(std::vector<Bond*>::iterator bt=bondVec.begin(); bt!=bondVec.end(); bt++){

			// check which pt1 or pt2 in (*bt) is the center, namely (*pt)
			bool is_pt1 = false;	// true when (*pt1) is the center
			if( this == (*bt)->getPt1() ){
				is_pt1 = true;
			}

		    	//Vec xi = (*bt)->getXi(is_pt1);
		    	//K += dyadicProduct(xi, xi)*(*bt)->getParticleVolume(is_pt1)*(*bt)->getWeight();
			K = K + (*bt)->getMicroK(is_pt1);

		} // end bond

		//// for numerical purpose, to be deleted later
		//K = 1.0/(horizonSize*horizonSize)*K;
		//
	 	//// inverse of matrix K
		//Kinv = K.getInvs()/(horizonSize*horizonSize);

		Kinv = inv(K);

	} // end calcParticleKinv()


	void PeriParticle::checkParticleAlive(REAL stretch_limit){

		int num_bonds = 0;	// the number of alive bonds
		for(std::vector<Bond*>::iterator bt=bondVec.begin(); bt!=bondVec.end(); bt++){

		    if( (*bt)->getIsAlive() ){	
			REAL bond_length = (*bt)->calcCurrentLength();

			REAL init_length = (*bt)->getInitLength();
			REAL stretch = ( bond_length - init_length )/init_length;
			
			if(stretch > stretch_limit || stretch < -2.0 ){
			    (*bt)->setAliveFalse();
			}
			else{
			    num_bonds++;
			}

		    } // if alive
	 	} // end bond

	 	// disable a particle
		if(num_bonds < 1){	// as rigid particle
		    isAlive = false;
			std::cout << "A particle is disabled due to the lack of bond" << std::endl;
		}

 	} // end checkParticleAlive()


	void PeriParticle::calcParticleStress(const PeriMaterialModel* periMaterial){


 		if( !isAlive ) {	// not alive
			sigma = zeros(3,3);}
	    else{
	    	// calculate deformation gradient tensor at current and half step
	    	Matrix N(3,3);	// matrix N at n+1 step
	    	Matrix N_half(3,3);	// matrix N at n+1/2 step
	    	Matrix N_deltaU(3,3); 
			// matrix N, corresponding to \mathbf{u}^{n+1} - \mathbf{u}^{n}, 
			// used to calculate \nabla (\mathbf{u}^{n+1} - \mathbf{u}^{n})

			for(std::vector<Bond*>::iterator bt=bondVec.begin(); bt!=bondVec.end(); bt++) {
				// check which pt1 or pt2 in (*bt) is the center, namely (*pt)
		        bool is_pt1 = false;	// true when (*pt1) is the center
		       	if( this == (*bt)->getPt1() ){
					is_pt1 = true;
		    	}

				bool bondIsAlive = (*bt)->getIsAlive();

				N = N + (*bt)->getMicroN(is_pt1,bondIsAlive);

				N_half = N_half + (*bt)->getMicroNHalf(is_pt1, bondIsAlive, TIMESTEP);
	
				N_deltaU = N_deltaU + (*bt)->getMicroNDeltaU(is_pt1, bondIsAlive, TIMESTEP);

				//if((*bt)->getIsAlive()){
	
				//	N += (*bt)->getMicroN(is_pt1);

				//	N_half += (*bt)->getMicroNHalf(is_pt1, TIMESTEP);
	
				//	N_deltaU += (*bt)->getMicroNDeltaU(is_pt1, TIMESTEP);
	
			
		  //  	}

			} // end bond

			deformationGradient = N*Kinv;
			deformationGradientHalf = N_half*Kinv;
			REAL eps = 1.0e-2;
			if(det(deformationGradient)<eps || det(deformationGradientHalf)<eps ){
				// calculate the determinant of deformationGraident and deformationGradientHalf, 
				// if the determinants are too small, then this particle is disabled, isAlive = false
				isAlive = false;	// disabled particle
				sigma = zeros(3,3);
				std::cout << "A particle is disabled because det[F] < 0.0" << std::endl;
	    	}
	    	else{
				if(periMaterial->getTypeConstitutive() == 1) {
					// Linear Elasticity, for testing purpose
					Matrix identity3x3(3,3);
					identity3x3(1,1) = 1; identity3x3(2,2) = 1; identity3x3(3,3) = 1;
					Matrix dudx = (deformationGradient - identity3x3)*(inv(deformationGradient));
					Matrix voight_strain(6,1);
					voight_strain(1,1) = dudx(1,1); 
					voight_strain(2,1) = dudx(2,2);
					voight_strain(3,1) = dudx(3,3);
					voight_strain(4,1) = dudx(2,3) + dudx(3,2);
					voight_strain(5,1) = dudx(1,3) + dudx(3,1); 
					voight_strain(6,1) = dudx(1,2) + dudx(2,1);
					Matrix voight_sigma = periMaterial->getTangentModulus()*voight_strain;
					sigma(1,1) = voight_sigma(1,1); sigma(2,2) = voight_sigma(2,1); 
					sigma(3,3) = voight_sigma(3,1); sigma(2,3) = voight_sigma(4,1);
					sigma(1,3) = voight_sigma(5,1); sigma(1,2) = voight_sigma(6,1); 
					sigma(2,1) = sigma(1,2); 
					sigma(3,1) = sigma(1,3);
					sigma(3,2) = sigma(2,3);
				}else if(periMaterial->getTypeConstitutive() == 2)
				{
	    			// calculate G, \nabla \Delta \bf{u}
					Matrix G = N_deltaU*Kinv*inv(deformationGradientHalf);
					Matrix Gsymm = 0.5*(G+trans(G));	// symmetric part of G
					Matrix Gskew = 0.5*(G-trans(G));	// skew part of G
	
					Matrix voight_Gsymm(6,1);
					voight_Gsymm(1,1) = Gsymm(1,1); voight_Gsymm(2,1) = Gsymm(2,2); 
					voight_Gsymm(3,1) = Gsymm(3,3); voight_Gsymm(4,1) = Gsymm(2,3);
					voight_Gsymm(5,1) = Gsymm(1,3); voight_Gsymm(6,1) = Gsymm(1,2);  
		
					Matrix voight_delta_sigma = periMaterial->getTangentModulus()*voight_Gsymm;
		
					Matrix delta_sigma(3,3);
					delta_sigma(1,1) = voight_delta_sigma(1,1); delta_sigma(2,2) = voight_delta_sigma(2,1); 
					delta_sigma(3,3) = voight_delta_sigma(3,1); delta_sigma(2,3) = voight_delta_sigma(4,1);
					delta_sigma(1,3) = voight_delta_sigma(5,1); delta_sigma(1,2) = voight_delta_sigma(6,1); 
					delta_sigma(2,1) = delta_sigma(1,2); 
					delta_sigma(3,1) = delta_sigma(1,3);
					delta_sigma(3,2) = delta_sigma(2,3);
		
					Matrix identity3x3(3,3);
					identity3x3(1,1) = 1; identity3x3(2,2) = 1; identity3x3(3,3) = 1;
					
					Matrix Q = identity3x3+inv(identity3x3-0.5*Gskew)*Gskew;
					Matrix trial_sigma = Q*sigma*trans(Q)+delta_sigma; 
		
					 // calculate deviatoric trial stress
					Matrix deviatoric_trial_sigma = trial_sigma;
					REAL trace_trial_sigma = trial_sigma(1,1)+trial_sigma(2,2)+trial_sigma(3,3);
					deviatoric_trial_sigma(1,1) = trial_sigma(1,1) - 1.0/3.0*trace_trial_sigma;
					deviatoric_trial_sigma(2,2) = trial_sigma(2,2) - 1.0/3.0*trace_trial_sigma;
					deviatoric_trial_sigma(3,3) = trial_sigma(3,3) - 1.0/3.0*trace_trial_sigma;
	
					REAL L2norm_deviatoric_trial_sigma = 0.0;
					L2norm_deviatoric_trial_sigma = deviatoric_trial_sigma(1,1)*deviatoric_trial_sigma(1,1) 
								      + deviatoric_trial_sigma(2,2)*deviatoric_trial_sigma(2,2)
								      + deviatoric_trial_sigma(3,3)*deviatoric_trial_sigma(3,3)
								      + 2.0*( deviatoric_trial_sigma(1,2)*deviatoric_trial_sigma(1,2) )
								      + 2.0*( deviatoric_trial_sigma(1,3)*deviatoric_trial_sigma(1,3) )
								      + 2.0*( deviatoric_trial_sigma(2,3)*deviatoric_trial_sigma(2,3) );
					L2norm_deviatoric_trial_sigma = sqrt(L2norm_deviatoric_trial_sigma);
	
					REAL Aphi = periMaterial->getAphi();
					REAL Bphi = periMaterial->getBphi();
					REAL cn = isv(1,1); //?
					REAL f_trial = L2norm_deviatoric_trial_sigma-(Aphi*cn-Bphi*1.0/3.0*trace_trial_sigma);
					
					if(f_trial < 0 ){	// elasticity
					    sigma = trial_sigma;
					    isv(1,1) = cn;
					}
					else{	// plasticity
	
					    REAL Bpsi = periMaterial->getBpsi();
					    REAL Hc = periMaterial->getHchi();
					    REAL KBulk = periMaterial->getKBulk();
					    REAL mu = periMaterial->getMu();
					    REAL delta_gamma = f_trial/( 2.0*mu + KBulk*Bphi*Bpsi + Hc*Aphi*Aphi);
					    sigma = trial_sigma - delta_gamma*( KBulk*Bpsi*identity3x3+2.0*mu*deviatoric_trial_sigma/L2norm_deviatoric_trial_sigma);
					    isv(1,1) = cn+delta_gamma*Hc*Aphi;
					}
				}

			}	// alive particle


		} //
	

  	} // end calcParticleStress()


	void PeriParticle::calcParticleAcceleration(const REAL density){

	    acceleration = 0.0;

	    Matrix acceleration_matrix(3,1);
	    Matrix xi_ik_matrix(3,1);
	    Matrix Pi;
	    Matrix Pk;
	    if(isAlive){
			for(std::vector<Bond*>::iterator bt=bondVec.begin(); bt!=bondVec.end(); bt++){

				PeriParticle* pti;
				PeriParticle* ptk;
				if( this == (*bt)->getPt1() ){
		    		pti = (*bt)->getPt1();
					ptk = (*bt)->getPt2();
				}
				else{
		    		pti = (*bt)->getPt2();
					ptk = (*bt)->getPt1();
				}

		
				// Piola Kirchoff stress of particle i
				Pi = det(pti->deformationGradient)*pti->sigma*inv( trans(pti->deformationGradient) );
				// Piola Kirchoff stress of particle k
				Pk = det(ptk->deformationGradient)*ptk->sigma*inv( trans(ptk->deformationGradient) );

				Vec xi_ik = ptk->initPosition - pti->initPosition;
				xi_ik_matrix(1,1) = xi_ik.getx(); 
				xi_ik_matrix(2,1) = xi_ik.gety();
				xi_ik_matrix(3,1) = xi_ik.getz();

				acceleration_matrix = acceleration_matrix + (*bt)->getWeight()*
					( Pi*(pti->Kinv) + Pk*(ptk->Kinv) )*xi_ik_matrix*ptk->particleVolume;


			} // end bond

			acceleration_matrix = acceleration_matrix/density;
			acceleration.setx(acceleration_matrix(1,1));
			acceleration.sety(acceleration_matrix(2,1));
			acceleration.setz(acceleration_matrix(3,1));

		}	// alive particle
	    

	} // end calcParticleAcceleration()


	void PeriParticle::updateDisplacement(const REAL dt){

	    velocityHalf = velocity + 0.5*acceleration*dt;
	    displacement += velocityHalf*dt;

	} // end updateDisplacement()


	void PeriParticle::updateVelocity(const REAL dt){

	    velocity = velocityHalf + 0.5*acceleration*dt;

	} // end updateVelocity()

	void PeriParticle::initial(){
	    displacement = 0.0;
	    velocity = 0.0;
		velocityHalf = 0.0;
	    acceleration = 0.0;
		sigma = zeros(3,3);
		deformationGradient = zeros(3,3);
		deformationGradientHalf = zeros(3,3);

	} // end initial()





} // end periDynamics

