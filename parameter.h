// #define marcos are placed in namespace dem nominally for convenience,
// they are not contained by any namespace, i.e., they are preprocessed.
 
#include "realtypes.h"
#include <fstream>

#ifndef PARAMETER_H
#define PARAMETER_H

namespace periDynamics { 

///////////////////////////////////////////////////////////////////////////////////////
// Part A: These parameters do not change frequently
///////////////////////////////////////////////////////////////////////////////////////

// PI value
extern const REAL PI;

// Gravitational acceleration
extern const REAL G;

// EPS (NOT float point relative precision, eps), problem domain dependent
extern const REAL EPS;

// relative overlap between particles
extern const REAL MINOVERLAP;
extern const REAL MAXOVERLAP;

// macro to toggle on/off MEPS
#define MEASURE_EPS

// measurable absolute overlap precision between particles
extern const REAL MEPS;

// random number seed
extern long idum;

// macro to toggle on/off random shape for each particle
//#define RANDOM_SHAPE

// particle material property
extern const REAL YOUNG;  
extern const REAL POISSON;      
extern const REAL Gs;   

// membrane particle material property
extern const REAL memYOUNG;
extern const REAL memPOISSON;

// other global variables
extern std::ofstream g_debuginf;
extern std::ofstream g_timeinf;
extern int g_iteration;

// output field width and precision
extern const int OWID;
extern const int OPREC;

///////////////////////////////////////////////////////////////////////////////////////
// Part B: These parameters may change frequently and can be easily edited in main.cpp
///////////////////////////////////////////////////////////////////////////////////////

// number of OpenMP threads
extern int  NUM_THREADS;

// 1. time integration method 
extern REAL TIMESTEP;
extern REAL MASS_SCL;
extern REAL MNT_SCL;
extern REAL GRVT_SCL;
extern REAL DMP_F;
extern REAL DMP_M;

// 2. normal damping and tangential friction
extern REAL DMP_CNT;
extern REAL FRICTION;
extern REAL BDRYFRIC;
extern REAL COHESION;

// 3. boundary displacement rate
extern REAL COMPRESS_RATE;
extern REAL RELEASE_RATE;
extern REAL PILE_RATE;
extern REAL STRESS_ERROR;

} // namespace periDynamics

#endif
