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
	char		tag[6]		= "petsc";
	int		nx[2]		= { 18, 18 };
	bool		vert[2]		= { false, true };
	bool		horiz[2]	= { true, false };
	bool		periodic[2]	= { false, false };
	double		min[2]		= { 0.0, 0.0 };
	double		max[2]		= { 1.0, 1.0 };
	Mesh*		vMesh		= new Mesh( "v-mesh", nx, "legendre", 5, min, max, periodic );
	Mesh*		hMesh		= new Mesh( "h-mesh", nx, "legendre", 3, min, max, periodic );
	BCs*		vbbcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*		v1bcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*		v2bcs		= new BCs( vert, vert, horiz, horiz, vMesh, 2*vMesh->nVertsTotal - v1bcs->size[0] - v1bcs->size[1] );
	BCs*		pbcs		= new BCs( false, false, false, false, vMesh, 0 );
	BCs*		hbcs		= new BCs( false, false, false, false, hMesh, 2*v2bcs->vecOffset );
	Field*		velTop		= new Field( "vel-top", vMesh, 2, v1bcs );
	Field*		velBot		= new Field( "vel-bot", vMesh, 2, v2bcs );
	Field*		velBar		= new Field( "vel-bar", vMesh, 2, vbbcs );
	Field*		etaInt		= new Field( "eta-int", hMesh, 1, hbcs );
	Field*		presSurf	= new Field( "pres-surf", vMesh, 1, pbcs );
	Field*		velTopPrev	= new Field( "vel-top-prev", vMesh, 2, v1bcs );
	Field*		velBotPrev	= new Field( "vel-bot-prev", vMesh, 2, v2bcs );
	Field*		velBarPrev	= new Field( "vel-bar-prev", vMesh, 2, vbbcs );
	Field*		etaIntPrev	= new Field( "eta-int-prev", hMesh, 1, hbcs );
	//double		dt		= 0.2/(3.0*nx[0]);
	double		dt		= 0.001;
	double		time		= 0.0;
	int		step		= atoi( argv[1] );
	int		dumpEvery	= 1;
	Params		topParams;
	Params		botParams;
	Params		barParams;
	Field*		vFields[4];
	Field*		hFields[1];
	BarotropicEqn*	bt;
	SW2L*		sw;
	Vec		G;

	PetscInitialize( &argc, &argv, (char)0, tag );

	topParams.f0    = botParams.f0    = barParams.f0    = 7.3e+1;
	topParams.beta  = botParams.beta  = barParams.beta  = 2.0e+1;
	topParams.nu    = botParams.nu    = barParams.nu    = 3.3e-4;
	topParams.gamma = 0.0;
	botParams.gamma                   = barParams.gamma = 0.0;
	topParams.pFac  = botParams.pFac  = barParams.pFac  = 1.0;
	topParams.hFac  = 0.0;
	botParams.hFac                    = barParams.hFac  = 2.0e+1;
	topParams.bFac  = 0.0;
	botParams.bFac                    = barParams.bFac  = 0.0;
	topParams.tau                     = barParams.tau   = 0.1;
	botParams.tau   = 0.0;
	topParams.k                       = barParams.k     = M_PI;
	botParams.k     = 0.0;
	topParams.H	= 1.05;
	botParams.H 	= 3.95;
	barParams.H     = 5.00;
	topParams.svv   = botParams.svv   = barParams.svv   = 9;

	vFields[0] = velTop;
	vFields[1] = velBot;
	vFields[2] = velBar;
	vFields[3] = presSurf;
	hFields[0] = etaInt;
	
	vMesh->Save();
	hMesh->Save();

	if( step > 0 ) {
		velTop->Read( step );
		velBot->Read( step );
		etaInt->Read( step );
		presSurf->Read( step );
	}

	/* system setup */
	bt = new BarotropicEqn( velBar, presSurf, dt, &barParams, NULL );
	sw = new SW2L( velTop, velBot, etaInt, velTopPrev, velBotPrev, etaIntPrev, presSurf, NULL, &topParams, &botParams, dt );

	time += dt;

	bt->Setup( 1 );
	sw->Assemble( 1 );

	/* solve first order */
	bt->AssembleGFirstOrder( &topParams, &botParams, velTop, velBot, etaInt, &G );
	bt->CalcVelBar( &topParams, &botParams, velTop, velBot, etaInt, velBar );
	bt->Solve( NULL, G );
	sw->Solve( 1 );
	VecDestroy( G );

	WriteXDMFHeader( step + 1 );
	WriteXDMF( vFields, 4, step + 1, time, dt );
	WriteXDMF( hFields, 1, step + 1, time, dt );
	WriteXDMFFooter( step + 1 );

	bt->Setup( 2 );
	sw->Assemble( 2 );
	for( int step_i = step + 2; step_i <= step + 40000; step_i++ ) {
		time += dt;
		cout << "solving for step: " << step_i << "\time: " << time << endl;
		/* solve second order */
		bt->AssembleGSecondOrder( &topParams, &botParams, velTop, velBot, etaInt, velTopPrev, velBotPrev, etaInt, &G );
		bt->CalcVelBar( &topParams, &botParams, velTop, velBot, etaInt, velBar );
		bt->CalcVelBar( &topParams, &botParams, velTopPrev, velBotPrev, etaIntPrev, velBarPrev );
		bt->Solve( velBarPrev, G );
		sw->Solve( 2 );
		VecDestroy( G );
		if( step_i%dumpEvery == 0 ) {
			WriteXDMFHeader( step_i );
			WriteXDMF( vFields, 4, step_i, time, dt );
			WriteXDMF( hFields, 1, step_i, time, dt );
			WriteXDMFFooter( step_i );
		}
	}

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
	delete velTopPrev;
	delete velBotPrev;
	delete velBarPrev;
	delete etaIntPrev;
	delete bt;
	delete sw;

	PetscFinalize();

	return EXIT_SUCCESS;
}
