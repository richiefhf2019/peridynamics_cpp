#include "PeriMaterialModel.h"
#include <math.h>


namespace periDynamics{

  	void PeriMaterialModel::calcMaterialParameters() {
	
		lambda=Poisson*Young/(1.0+Poisson)/(1.0-2.0*Poisson);
   		mu=Young/2.0/(1.0+Poisson);
       	kBulk=Young/(3.0*(1.0-2.0*Poisson));

       	tangentModulus(1,1) = lambda+2.00*mu;
       	tangentModulus(1,2) = lambda;
       	tangentModulus(1,3) = lambda;
       	tangentModulus(2,1) = lambda; 
       	tangentModulus(2,2) = lambda+2.00*mu ;
       	tangentModulus(2,3) = lambda; 
       	tangentModulus(3,1) = lambda; 
       	tangentModulus(3,2) = lambda; 
       	tangentModulus(3,3) = lambda+2.00*mu;
       	tangentModulus(4,4) = mu; 
       	tangentModulus(5,5) = mu; 
      	tangentModulus(6,6) = mu;
      
      	hchi = 0.00;
     	chi=208848.10;	//kPa
      	c=208848.10;	//kPa
      	phi=0.7386630; //radians
      	psi=0.7386630; //radians

      	kappa=-5.e4; 	//kPa
      	rEllip=1.0;   	//dimensionless

      	beta=-1.0;   	//TC
		bondStretchLimit = 0.45;

		Aphi = 2*sqrt(6.0)*cos(phi)/(3.0+beta*sin(phi));
		Bphi = 2*sqrt(6.0)*sin(phi)/(3.0+beta*sin(phi));

		Apsi = 2*sqrt(6.0)*cos(psi)/(3.0+beta*sin(psi));
		Bpsi = 2*sqrt(6.0)*sin(psi)/(3.0+beta*sin(psi));

    } // calcMaterialParameters()


} // end periDynamics
