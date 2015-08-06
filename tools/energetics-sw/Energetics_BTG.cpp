#include <cmath>
#include <iostream>
#include <fstream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SW.h"

using namespace std;
using std::string;

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 18, 36 };
	bool			periodic[2]	= { false, false };
	bool			noSlip[2]	= { true, true };
	double			min[2]		= { 0.0, 0.0 };
	double			max[2]		= { 1.0, 2.0 };
	Mesh*			vMesh		= new Mesh( "vMesh", nx, "legendre", 5, min, max, periodic );
	Mesh*			pMesh		= new Mesh( "pMesh", nx, "legendre", 3, min, max, periodic );
	BCs*			vbcs		= new BCs( noSlip, noSlip, noSlip, noSlip, vMesh, 0 );
	BCs*			pbcs		= new BCs( false, false, false, false, pMesh, 2*vMesh->nVertsTotal - vbcs->size[0] - vbcs->size[1] );
	Field*			vel		= new Field( "vel", vMesh, 2, vbcs );
	Field*			phi		= new Field( "phi", pMesh, 1, pbcs );
	double			dt		= 0.5/(3*nx[0]);
	int			skip		= atoi( argv[1] );
	int			end		= atoi( argv[2] );
	ofstream		file;
	SWParams		params;
	ShallowWaterEqn* 	sw;
	Energetics_SW*		en;

	PetscInitialize( &argc, &argv, (char)0, tag );

	params.L     = 1.0;
	params.H     = 1.0;
	params.U     = 1.0;
	params.rho   = 1.0;
	params.f0    = 50.0;
	params.beta  = 20.0;
	params.g     = 15.0;
	params.gamma = 0.05;
	params.tau   = 0.19569471624266144*0.8;
	params.kws   = M_PI;
	params.nu    = 3.0e-4;

	sw = new ShallowWaterEqn( vel, phi, dt, 4, NULL, &params );
	en = new Energetics_SW( sw );

	file.open( "sw.en" );

	for( int step_i = skip; step_i <= end; step_i += skip ) {
		file << "\t" << step_i;
		file << "\t" << step_i*dt;
		file << "\t" << en->CalcKE( step_i );
		file << "\t" << en->CalcPE( step_i );
		file << "\t" << en->CalcPower_WindStress( step_i );
		file << "\t" << en->CalcPower_Friction( step_i );
		file << "\n";
	}

	file.close();

	delete sw;
	delete en;
	delete vMesh;
	delete pMesh;
	delete vbcs;
	delete pbcs;
	delete vel;
	delete phi;

	PetscFinalize();

	return EXIT_SUCCESS;
}
