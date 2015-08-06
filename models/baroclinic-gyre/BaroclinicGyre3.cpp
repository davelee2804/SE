#include <iostream>
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
two layer baroclinic gyre model

reference:
	J. J. Nauw and H. A. Dijkstra (2001) "The origin of low-frequency variability 
	of double-gyre wind-driven flows" J. Mar. Res, 59, 567-597

	params:
		U	1.0		m/s
		L	1.0x10^6	m
		H	1.0x10^3	m
		H_1	7.0x10^2	m
		H_2	3.3x10^3	m
		f_0	1.0x10^-4	1/s
		beta	1.8x10^-11	1/ms
		rho_0	1.0x10^3	kg/m^3
		rho_1	1.000x10^3	kg/m^3
		rho_2	1.002x10^3	kg/m^3
		tau	0.1		kg/ms^2
		nu	3x102		m^2/s
		gamma	5.0x10^-8	1/s
*/

int main( int argc, char** argv ) {
	char		tag[6]		= "petsc";
	int		nx[2]		= { 8, 16 };
	bool		noSlip[2]	= { true, true };
	bool		periodic[2]	= { false, false };
	double		min[2]		= { 0.0, 0.0 };
	double		max[2]		= { 1.0, 2.0 };
	Mesh*		vMesh		= new Mesh( "v-mesh", nx, "legendre", 8, min, max, periodic );
	Mesh*		hMesh		= new Mesh( "h-mesh", nx, "legendre", 6, min, max, periodic );
	BCs*		vbcs		= new BCs( noSlip, noSlip, noSlip, noSlip, vMesh, 0 );
	BCs*		pbcs		= new BCs( false, false, false, false, vMesh, 0 );
	BCs*		hbcs		= new BCs( false, false, false, false, hMesh, 0 );
	Field*		velTop		= new Field( "vel-top", vMesh, 2, vbcs );
	Field*		velBot		= new Field( "vel-bot", vMesh, 2, vbcs );
	Field*		velBar		= new Field( "vel-bar", vMesh, 2, vbcs );
	Field*		etaInt		= new Field( "eta-int", hMesh, 1, hbcs );
	Field*		presSurf	= new Field( "pres-surf", vMesh, 1, pbcs );
	Field*		velTopPrev	= new Field( "vel-top-prev", vMesh, 2, vbcs );
	Field*		velBotPrev	= new Field( "vel-bot-prev", vMesh, 2, vbcs );
	Field*		velBarPrev	= new Field( "vel-bar-prev", vMesh, 2, vbcs );
	Field*		etaIntPrev	= new Field( "eta-int-prev", hMesh, 1, hbcs );
	Field*		topography	= new Field( "bot-topo", hMesh, 1, NULL );
	double		dt		= 0.25*(max[0] - min[0])/(nx[0]*8);
	double		time		= 0.0;
	Params		topParams;
	Params		botParams;
	Field*		vFields[4];
	Field*		hFields[1];
	TwoLayerEqn*	eqn;
	MomentumEqn*	tl;

	PetscInitialize( &argc, &argv, (char)0, tag );

	topParams.f0    = botParams.f0    = 1.0e+2;
	topParams.beta  = botParams.beta  = 1.8e+1;
	topParams.nu    = botParams.nu    = 3.0e-4;
	topParams.gamma = 0.0;
	botParams.gamma = 5.0e-2;
	topParams.pFac  = botParams.pFac  = 1.96e+1;
	topParams.hFac  = 0.0;
	botParams.hFac  = 1.96e+1;
	topParams.bFac  = 0.0;
	botParams.bFac  = 1.96e+1;
	topParams.tau   = 0.1;	/* 0.05?? */
	botParams.tau   = 0.0;
	topParams.k     = M_PI;
	botParams.k     = 0.0;
	topParams.H	= 0.7;
	botParams.H 	= 3.3;
	topParams.svv   = botParams.svv    = vMesh->el->N/2;

	vFields[0] = velTop;
	vFields[1] = velBot;
	vFields[2] = velBar;
	vFields[3] = presSurf;
	hFields[0] = etaInt;
	vMesh->Save();
	hMesh->Save();

	/* initial guess for the top layer */
	tl = new MomentumEqn( velTop, etaInt, velTopPrev, etaIntPrev, NULL, NULL, &topParams, dt );
	tl->Assemble( 1 );
	tl->Solve( 1 );
	tl->Update();
	delete tl;
	WriteXDMFHeader( 0 );
	WriteXDMF( &velTop, 1, 0, time, dt );
	WriteXDMFFooter( 0 );

	eqn = new TwoLayerEqn( velTop, velBot, etaInt, velBar, presSurf, topography, velTopPrev, velBotPrev, etaIntPrev, velBarPrev, &topParams, &botParams, dt );
	time += dt;
	eqn->Setup( 1 );
	eqn->Solve( 1 );
	
	WriteXDMFHeader( 1 );
	WriteXDMF( vFields, 4, 1, time, dt );
	WriteXDMF( hFields, 1, 1, time, dt );
	WriteXDMFFooter( 1 );

	eqn->Setup( 2 );
	for( int step_i = 2; step_i <= 100; step_i++ ) {
		time += dt;
		eqn->Solve( 2 );
		WriteXDMFHeader( step_i );
		WriteXDMF( vFields, 4, step_i, time, dt );
		WriteXDMF( hFields, 1, step_i, time, dt );
		WriteXDMFFooter( step_i );
	}

	delete vMesh;
	delete hMesh;
	delete vbcs;
	delete pbcs;
	delete hbcs;
	delete velTop;
	delete velBot;
	delete velBar;
	delete etaInt;
	delete presSurf;
	delete velTopPrev;
	delete velBotPrev;
	delete velBarPrev;
	delete etaIntPrev;
	delete topography;
	delete eqn;

	PetscFinalize();

	return EXIT_SUCCESS;
}
