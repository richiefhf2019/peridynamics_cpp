#pragma once
#ifndef GLOBFUNCS_H
#define GLOBFUNCS_H

#include "Vec.h"
#include "Matrix.h"
#include "realtypes.h"
#include <iostream>

namespace periDynamics{

	void gauss3D(const int nip, Matrix& gp_loc3D, Matrix& gp_weigth3D);
	
	void gauss1D(const int nintElem, Matrix& gp_loc1D, Matrix& gp_weight1D );
	
	void invert(const Matrix& a, Matrix& a_inv, REAL& determinant);

	void shp3d(const REAL, const REAL, const REAL, Matrix&, Matrix&, REAL&);
	
	REAL windowfunction(REAL length);

	Matrix dyadicProduct(const Vec&, const Vec&);

} // end periDynamics

#endif
