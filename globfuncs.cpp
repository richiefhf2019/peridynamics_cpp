#include "globfuncs.h"
#include <stdlib.h>

namespace periDynamics{


    //***********************************************************
    //
    //   Use 1-D Gauss data to generate quadrature data for 
    //   a 3-d brick element (cube)
    //
    //** input parameter:
    //   nip:   the number of Gauss point in one direction
    //   nint:  the total number of Gauss points per element
    //          nintElem = nip * nip * nip              
    //
    //   (Note: nint is the name of Math intrinsic function )
    //
    //   gp_loc3D: the local coordinate of the Gauss quadrature;
    //   gp_weight3D: the weight of the Gauss quadrature;
    //
    //   Shaofan Li,  August, 1998
    //
    //*************************************************************
    void gauss3D( const int nip, Matrix& gp_loc3D, Matrix& gp_weight3D ) {
		// gp_loc3D is already 3x8, and gp_weight3D is already 1x8

		gp_loc3D = zeros(3,8);
		gp_weight3D = zeros(1,8);

		Matrix gp_loc1D;
		Matrix gp_weight1D;

		gauss1D( nip, gp_loc1D, gp_weight1D);

		int ncount = 0;
		for(int ig = 0; ig < nip; ig++){
			for(int jg = 0; jg < nip; jg++){
				for(int kg = 0; kg < nip; kg++){
					ncount++;
					gp_weight3D(1,ncount) = gp_weight1D(1,ig+1)*gp_weight1D(1,jg+1)*gp_weight1D(1,kg+1);
					gp_loc3D(1,ncount) = gp_loc1D(1,ig+1);
					gp_loc3D(2,ncount) = gp_loc1D(1,jg+1);
					gp_loc3D(3,ncount) = gp_loc1D(1,kg+1);
				}
			}
		}

    } // end gauss3D()

   void gauss1D(const int nintElem, Matrix& s, Matrix& w ) {
		// s and w now are not given spaces
		// subroutine to give gaussian pts (up to 10) of 1D
		// for intergration over -1 to 1 !!!!!!

		if( nintElem > 10 ){
			std::cout << "nintElem > 10 in subroutine of gaussian. STOP!" << std::endl;
			std::cout << "nintElem = " << nintElem << std::endl;
			exit(1);
			}

		s = zeros(1,nintElem);
		w = zeros(1,nintElem);
		switch (nintElem) {
			case 1:
				s(1,1) = 0.0;
				w(1,1) = 2.0;
	 			break;
			case 2:
	 			s(1,1) = -0.5773502691896260;
      			s(1,2) = -s(1,1);
      			w(1,1) = 1.0;
      			w(1,2) = 1.0;
				break;
			case 3:
				s(1,1) = -0.7745966692414830;
      			s(1,2) = 0.0;
      			s(1,3) = -s(1,1);
      			w(1,1) = 0.5555555555555560;
      			w(1,2) = 0.8888888888888890;
      			w(1,3) = w(1,1);
			break;
			case 4:
				s(1,1) = -0.8611363115940530;
      			s(1,2) = -0.3399810435848560;
      			s(1,3) = -s(1,2);
      			s(1,4) = -s(1,1);
      			w(1,1) = 0.3478548451374540;
      			w(1,2) = 0.6521451548625460;
      			w(1,3) = w(1,2);
      			w(1,4) = w(1,1);
			break;
			case 5:
				s(1,1) = -0.9061798459386640;
      			s(1,2) = -0.5384693101056830;
      			s(1,3) = 0.;
      			s(1,4) = -s(1,2);
      			s(1,5) = -s(1,1);
     		 	w(1,1) = 0.2369368850561890;
      			w(1,2) = 0.4786386704993660;
      			w(1,3) = 0.5688888888888890;
      			w(1,4) = w(1,2);
      			w(1,5) = w(1,1);
			break;
			default:
			std::cout << "nintElem larger than 5 is not defined..." << std::cout;
			break;

	} // end switch

    } // end gauss1D()

	void shp3d(const REAL xi, const REAL eta, const REAL zeta, Matrix& xl, Matrix& shp, REAL& xsj) {

		REAL xim, xip, etam, etap, zetam, zetap;
		const REAL pt125 = 0.125;
		shp = zeros(4,8);
		xim = 1.0 - xi, etam = 1.0 - eta, zetam = 1.0 - zeta;
		xip = 1.0 + xi, etap = 1.0 + eta, zetap = 1.0 + zeta;

		// shape function (N^[1-8]) evaulated at point (xi, eta, zeta)
		shp(4,1) = pt125*xim*etam*zetam, shp(4,2) = pt125*xip*etam*zetam;
		shp(4,3) = pt125*xip*etap*zetam, shp(4,4) = pt125*xim*etap*zetam;
		shp(4,5) = pt125*xim*etam*zetap, shp(4,6) = pt125*xip*etam*zetap;
		shp(4,7) = pt125*xip*etap*zetap, shp(4,8) = pt125*xim*etap*zetap;

		// natural derivatives of shape functions evaluated at (xi, eta, zeta)
		shp(1,1)= -pt125*etam*zetam, shp(1,2)=  pt125*etam*zetam;
		shp(1,3)=  pt125*etap*zetam, shp(1,4)= -pt125*etap*zetam;
		shp(1,5)= -pt125*etam*zetap, shp(1,6)=  pt125*etam*zetap;
		shp(1,7)=  pt125*etap*zetap, shp(1,8)= -pt125*etap*zetap;
		
		shp(2,1)= -pt125*xim*zetam, shp(2,2)= -pt125*xip*zetam;
		shp(2,3)=  pt125*xip*zetam, shp(2,4)=  pt125*xim*zetam;
		shp(2,5)= -pt125*xim*zetap, shp(2,6)= -pt125*xip*zetap;
		shp(2,7)=  pt125*xip*zetap, shp(2,8)=  pt125*xim*zetap;
		 
		shp(3,1)= -pt125*xim*etam, shp(3,2)= -pt125*xip*etam;
		shp(3,3)= -pt125*xip*etap, shp(3,4)= -pt125*xim*etap;
		shp(3,5)=  pt125*xim*etam, shp(3,6)=  pt125*xip*etam;
		shp(3,7)=  pt125*xip*etap, shp(3,8)=  pt125*xim*etap;

        // loop to the find the jacobian matrix
        Matrix jac(3,3);
		for(int i = 1; i <= 3; i++)
			for(int j = 1; j <= 3; j++)
				for(int node = 1; node <= 8; node++) {
					jac(i,j) += shp(j,node)*xl(i,node);
				}
		
		xsj = det(jac);
	} // end shp3d


	Matrix dyadicProduct(const Vec& a, const Vec& b) {
	
		Matrix c = zeros(3,3);
		c(1,1) = a.getx()*b.getx();
		c(1,2) = a.getx()*b.gety();
		c(1,3) = a.getx()*b.getz();

		c(2,1) = a.gety()*b.getx();
		c(2,2) = a.gety()*b.gety();
		c(2,3) = a.gety()*b.getz();

		c(3,1) = a.getz()*b.getx();
		c(3,2) = a.getz()*b.gety();
		c(3,3) = a.getz()*b.getz();

		return c;

	}

} // end periDynamics
