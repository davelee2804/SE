#include <iostream>
#include <fstream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SW.h"
#include "2LRL.h"

using namespace std;
using std::string;

/*
two layer baroclinic gyre model - energetics test

reference:
	Holland, W. R. and Lin, L. B. (1975) "On the Generation of Mesoscale Eddies and their 
	Contribution to the Oceanic General Circulation. I. A Preliminary Numerical Experiment"
	J. Phys. Oceanogr. 5, 642-657

parameters:
	f_0	= 7.5e-5	s^-1
	beta	= 2.0e-11	m^-1s^-1
	nu	= 3.3e+2	m^2s^-1
	L	= 1.0e+6	m
	H	= 1.0e+3	m		( or 5.0e+3 ?? )
	g'	= 0.02		m.s^-2
	U	= 1		m.s^-1
	tau	= 0.1		kg.m^-1s^-2	( or 0.02 ?? )
	rho	= 1.0e+3	kg.m^-3
*/

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 32, 32 };
	bool			vert[2]		= { false, true };
	bool			horiz[2]	= { true, false };
	bool			periodic[2]	= { false, false };
	double			min[2]		= { -0.5, -0.5 };
	double			max[2]		= { +0.5, +0.5 };
	Mesh*			vMesh		= new Mesh( "v-mesh", nx, "legendre", 5, min, max, periodic );
	Mesh*			hMesh		= new Mesh( "h-mesh", nx, "legendre", 3, min, max, periodic );
	BCs*			vbbcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*			v1bcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*			v2bcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*			pbcs		= new BCs( false, false, false, false, vMesh, 0 );
	BCs*			hbcs		= new BCs( false, false, false, false, hMesh, 0 );
	Field*			velTop		= new Field( "vel-top", vMesh, 2, v1bcs );
	Field*			velBot		= new Field( "vel-bot", vMesh, 2, v2bcs );
	Field*			velBar		= new Field( "vel-bar", vMesh, 2, vbbcs );
	Field*			etaInt		= new Field( "eta-int", hMesh, 1, hbcs );
	Field*			presSurf	= new Field( "pres-surf", vMesh, 1, pbcs );
	double			dt		= 0.0005;
	int			skip		= atoi( argv[1] );
	int			end		= atoi( argv[2] );
	TwoLayerParams		params;
	TwoLayerSW*		sw;
	TwoLayerSW_Energetics*	en;
	ofstream		file;

	PetscInitialize( &argc, &argv, (char)0, tag );

	params.f     = 7.3e+1;
	params.beta  = 2.0e+1;
	params.nu    = 3.3e-4;
	//params.gamma = 0.0;
	params.p     = 1.0;
	params.g     = 2.0e+1;
	params.tau   = 0.1;
	params.k     = M_PI;
	params.H1    = 1.05;
	params.H2    = 3.95;

	/* system setup */
	sw = new TwoLayerSW( velTop, velBot, etaInt, presSurf, &params, dt );
	en = new TwoLayerSW_Energetics( sw );

	file.open( "twolayer_sw.en" );

	for( int step_i = skip; step_i <= end; step_i += skip ) {
		cout << "calculating energetics for time step: " << step_i << endl;
		file << "\t" << step_i;
		file << "\t" << step_i*dt;
		file << "\t" << en->CalcKE( step_i, 1 );
		file << "\t" << en->CalcKE( step_i, 2 );
		file << "\t" << en->CalcPE( step_i );
		file << "\t" << en->CalcPower_SurfPres( step_i, 1 );
		file << "\t" << en->CalcPower_SurfPres( step_i, 2 );
		file << "\t" << en->CalcPower_WindStress( step_i );
		//file << "\t" << en->CalcPower_Friction( step_i );
		file << "\t" << en->CalcPower_Viscosity( step_i, 1 );
		file << "\t" << en->CalcPower_Viscosity( step_i, 2 );
		file << "\n";
	}

	file.close();

	delete vMesh;
	delete hMesh;
	delete vbbcs;
	delete v1bcs;
	delete v2bcs;
	delete pbcs;
	delete hbcs;
	delete velTop;
	delete velBot;
	delete velBar;
	delete etaInt;
	delete presSurf;
	delete sw;
	delete en;

	PetscFinalize();

	return EXIT_SUCCESS;
}
