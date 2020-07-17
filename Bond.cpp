// Function Definitions

#include "Bond.h"
#include "PeriParticle.h"
namespace periDynamics{

	//-------------------------------------------------------------------------
	// Default Constructor
	Bond::Bond() {
		isAlive = true;
		weight = 0.0;
		initLength = 0.0;
		pt1 = NULL;
		pt2 = NULL;
	}
	
	// Overload Constructor
	Bond::Bond(REAL tmp_length, PeriParticle* tmp_pt1, PeriParticle* tmp_pt2):initLength(tmp_length), pt1(tmp_pt1), pt2(tmp_pt2){
	    isAlive = true;
     	} // Bond()
	
	// Destructor
	Bond::~Bond() {

	}

	//-------------------------------------------------------------------------
	// Accessor Functions
	bool Bond::getIsAlive() const {
		return isAlive;
	}

	REAL Bond::getWeight() const {
		return weight;
	}

	REAL Bond::getInitLength() const {
		return initLength;
	}

	REAL Bond::getParticleVolume(bool is_pt1) const {
	    if(is_pt1){
		return pt1->getParticleVolume();
	    }
	    else{
		return pt2->getParticleVolume();
	    }

	}

	Vec Bond::getXi(bool is_pt1) const {
		if(is_pt1){
		    return (pt2->getInitPosition()-pt1->getInitPosition());		
		}
		else{
		    return (pt1->getInitPosition()-pt2->getInitPosition());
		}
	} // end getXi()

	Vec Bond::getEta(bool is_pt1) const {
		if(is_pt1){
		    return (pt2->getInitPosition()+pt2->getDisplacement() - pt1->getInitPosition()-pt1->getDisplacement());		
		}
		else{
		    return (pt1->getInitPosition()+pt1->getDisplacement() - pt2->getInitPosition()-pt2->getDisplacement());
		}
	} // end getEta()

	Vec Bond::getEtaHalf(bool is_pt1, const REAL dt) const {
		if(is_pt1){
		    return (pt2->getInitPosition()+pt2->getDisplacement()-0.5*dt*pt2->getVelocityHalf()
				 - (pt1->getInitPosition()+pt1->getDisplacement()-0.5*dt*pt1->getVelocityHalf()) );		
		}
		else{
		    return (pt1->getInitPosition()+pt1->getDisplacement()-0.5*dt*pt1->getVelocityHalf() 
				 - (pt2->getInitPosition()+pt2->getDisplacement()-0.5*dt*pt2->getVelocityHalf()) );
		}
	} // end getEtaHalf()	

	Matrix Bond::getMicroK(const bool is_pt1) const {

		Vec xi;
		REAL volume;
		if(is_pt1){
		    xi = (pt2->getInitPosition()-pt1->getInitPosition());		
		    volume = pt1->getParticleVolume();
		}
		else{
		    xi = (pt1->getInitPosition()-pt2->getInitPosition());
		    volume = pt2->getParticleVolume(); 
		}

		return (dyadicProduct(xi, xi)*volume*weight);

	} // end getMicroK()

	Matrix Bond::getMicroN(const bool is_pt1, const bool bondIsAlive) const {

		Vec xi;
		Vec eta;
		REAL volume;
		if(bondIsAlive) {
			if(is_pt1){
			    xi = (pt2->getInitPosition()-pt1->getInitPosition());	
			    eta = (pt2->getInitPosition()+pt2->getDisplacement() - pt1->getInitPosition()-pt1->getDisplacement());	
			    volume = pt1->getParticleVolume();
			}
			else{
			    xi = (pt1->getInitPosition()-pt2->getInitPosition());
			    eta = (pt1->getInitPosition()+pt1->getDisplacement() - pt2->getInitPosition()-pt2->getDisplacement());
			    volume = pt2->getParticleVolume(); 
			}
		}
		else {
			if(is_pt1){
			    xi = (pt2->getInitPosition()-pt1->getInitPosition());	
			    //eta = (pt2->getInitPosition()+pt1->getDisplacement() - pt1->getInitPosition()-pt1->getDisplacement());
				eta = xi;
			    volume = pt1->getParticleVolume();
			}
			else{
			    xi = (pt1->getInitPosition()-pt2->getInitPosition());
			    //eta = (pt1->getInitPosition()+pt2->getDisplacement() - pt2->getInitPosition()-pt2->getDisplacement());
				eta = xi;
			    volume = pt2->getParticleVolume(); 
			}
		}

		return (dyadicProduct(eta, xi)*volume*weight);

	} // end getMicroN()

	Matrix Bond::getMicroNHalf(const bool is_pt1, const bool isBondAlive, const REAL dt) const {

		Vec xi;
		Vec eta;
		REAL volume;
		if(isBondAlive) {
			if(is_pt1){
			    xi = (pt2->getInitPosition()-pt1->getInitPosition());	
			    eta = (pt2->getInitPosition()+pt2->getDisplacement()-0.5*dt*pt2->getVelocityHalf()
					 - (pt1->getInitPosition()+pt1->getDisplacement()-0.5*dt*pt1->getVelocityHalf()) );	
			    volume = pt1->getParticleVolume();
			}
			else{
			    xi = (pt1->getInitPosition()-pt2->getInitPosition());
			    eta = (pt1->getInitPosition()+pt1->getDisplacement()-0.5*dt*pt1->getVelocityHalf() 
					 - (pt2->getInitPosition()+pt2->getDisplacement()-0.5*dt*pt2->getVelocityHalf()) );
			    volume = pt2->getParticleVolume(); 
			}
		}
		else {
			if(is_pt1){
			    xi = (pt2->getInitPosition()-pt1->getInitPosition());	
			    //eta = (pt2->getInitPosition()+pt2->getDisplacement()-0.5*dt*pt2->getVelocityHalf()
				//	 - (pt1->getInitPosition()+pt1->getDisplacement()-0.5*dt*pt1->getVelocityHalf()) );	
				eta = xi;
			    volume = pt1->getParticleVolume();
			}
			else{
			    xi = (pt1->getInitPosition()-pt2->getInitPosition());
			    //eta = (pt1->getInitPosition()+pt1->getDisplacement()-0.5*dt*pt1->getVelocityHalf() 
				//	 - (pt2->getInitPosition()+pt2->getDisplacement()-0.5*dt*pt2->getVelocityHalf()) );
				eta = xi;
			    volume = pt2->getParticleVolume(); 
			}

		}

		return (dyadicProduct(eta, xi)*volume*weight);

	} // end getMicroNHalf()

	Matrix Bond::getMicroNDeltaU(const bool is_pt1, const bool isBondAlive, const REAL dt) const {

		Vec xi;
		Vec eta;
		REAL volume;
		if(isBondAlive) {
			if(is_pt1){
			    xi = (pt2->getInitPosition()-pt1->getInitPosition());	
			    eta = dt*(pt2->getVelocityHalf() - pt1->getVelocityHalf());	
			    volume = pt1->getParticleVolume();
			}
			else{
			    xi = (pt1->getInitPosition()-pt2->getInitPosition());
			    eta = dt*(pt1->getVelocityHalf() - pt2->getVelocityHalf());	
			    volume = pt2->getParticleVolume(); 
			}
		}
		//commented out, because eta = 0.0 for this case, no contribution needs to be added
		//else {
		//	if(is_pt1){
		//	    xi = (pt2->getInitPosition()-pt1->getInitPosition());	
		//	    // eta = dt*(pt2->getVelocityHalf() - pt1->getVelocityHalf());	
		//		eta = 0.0;
		//	    volume = pt1->getParticleVolume();
		//	}
		//	else{
		//	    xi = (pt1->getInitPosition()-pt2->getInitPosition());
		//	    // eta = dt*(pt1->getVelocityHalf() - pt2->getVelocityHalf());	
		//		eta = 0.0;
		//	    volume = pt2->getParticleVolume(); 
		//	}
		//}

		return (dyadicProduct(eta, xi)*volume*weight);

	} // end getMicroNDeltaU()


	
	//-------------------------------------------------------------------------
	// Mutator Functions
	void Bond::setIsAlive(bool newisAlive) { 
		isAlive = newisAlive;
	}

	void Bond::setWeight(REAL newweight) {
		weight = newweight;
	}

	void Bond::setInitLength(REAL newinitLength) {
		initLength = newinitLength;
	}

	//-------------------------------------------------------------------------
	// Utility Functions

	REAL Bond::calcCurrentLength() {
		return vfabs( pt1->getInitPosition() + pt1->getDisplacement()
			        - pt2->getInitPosition() - pt2->getDisplacement() );
	}


} // end periDynamics
