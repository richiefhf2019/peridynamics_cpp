#include <iostream>
#include <string>
#include <vector>

// user defined head files
#include "PeriDomain.h"

using namespace std;
using namespace periDynamics;

// Main Function
int main(int argc, char * argv[]) {


	PeriDomain periCase;
	periCase.initial("blockTensionSmall.in");
	periCase.solve("blockTensionSmall.dat");
	//periCase.initial("linear.in");
	//periCase.solve("linear.dat");
	//periCase.readData("peridp.test");
	//periCase.calcParticleVolume();
	//periCase.calcHorizonSize();
	//periCase.constructNeighbor();
	//periCase.calcParticleKinv();
	//periCase.writeMeshCheckVolume("dp.dat");

	return 0;
}
