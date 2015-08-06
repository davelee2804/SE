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
	double		min[2]		= { -0.5, -0.5 };
	double		max[2]		= { +0.5, +0.5 };
	Mesh*		vMesh		= new Mesh( "v-mesh", nx, "legendre", 5, min, max, periodic );
	Mesh*		hMesh		= new Mesh( "h-mesh", nx, "legendre", 3, min, max, periodic );
	BCs*		v1bcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*		v2bcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*		pbcs		= new BCs( false, false, false, false, vMesh, 0 );
	BCs*		hbcs		= new BCs( false, false, false, false, hMesh, 0 );
	Field*		velTop		= new Field( "vel-top", vMesh, 2, v1bcs );
	Field*		velBot		= new Field( "vel-bot", vMesh, 2, v2bcs );
	Field*		etaInt		= new Field( "eta-int", hMesh, 1, hbcs );
	Field*		presSurf	= new Field( "pres-surf", vMesh, 1, pbcs );
	Field*		velTopPrev	= new Field( "vel-top-prev", vMesh, 2, v1bcs );
	Field*		velBotPrev	= new Field( "vel-bot-prev", vMesh, 2, v2bcs );
	Field*		etaIntPrev	= new Field( "eta-int-prev", hMesh, 1, hbcs );
	//double		dt		= 0.00025;
	double		dt		= 0.001;
	double		time		= 0.0;
	int		step		= atoi( argv[1] );
	int		dumpEvery	= 40;
	TwoLayerParams	params;
	Field*		vFields[3];
	Field*		hFields[1];
	TwoLayerSW*	sw;
        BarotropicEqn*  bt;
        BCs*            vbbcs           = new BCs( vert, vert, horiz, horiz, vMesh, 0 );
        Field*          velBar          = new Field( "vel-bar", vMesh, 2, vbbcs );
        Field*          velBarPrev      = new Field( "vel-bar-prev", vMesh, 2, vbbcs );
        Params          top;
        Params          bot;
        Params          bar;
        Vec             G;

	PetscInitialize( &argc, &argv, (char)0, tag );

	params.f     	= 7.3e+1;
	params.beta  	= 2.0e+1;
	params.nu    	= 3.3e-4;
	params.p 	= 1.0;
	params.g	= 2.0e+1;
	params.tau 	= 0.1;
	params.k 	= M_PI;
	params.H1	= 1.05;
	params.H2 	= 3.95;
        bar.f0 = 7.3e+1; bar.beta = 2.0e+1; bar.nu = 3.3e-4; bar.pFac = 1.0; bar.hFac = 2.0e+1; bar.bFac = 0.0; bar.tau = 0.1; bar.k = M_PI; bar.H = 5.00; bar.svv = 9;
        top.f0 = 7.3e+1; top.beta = 2.0e+1; top.nu = 3.3e-4; top.pFac = 1.0; top.hFac = 0.0000; top.bFac = 0.0; top.tau = 0.1; top.k = M_PI; top.H = 1.05; top.svv = 9;
        bot.f0 = 7.3e+1; bot.beta = 2.0e+1; bot.nu = 3.3e-4; bot.pFac = 1.0; bot.hFac = 2.03+1; bot.bFac = 0.0; bot.tau = 0.0; bot.k = 0.00; bot.H = 3.95; bot.svv = 9;

	vFields[0] = velTop;
	vFields[1] = velBot;
	vFields[2] = presSurf;
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
	sw = new TwoLayerSW( velTop, velBot, etaInt, presSurf, &params, dt );
	//sw->enull = true;
	//sw->pnull = true;
        bt = new BarotropicEqn( velBar, presSurf, dt, &bar, NULL );

	time += dt;

	/* solve first order */
        bt->Setup(1);
        bt->AssembleGFirstOrder( &top, &bot, velTop, velBot, etaInt, &G );
        bt->CalcVelBar( &top, &bot, velTop, velBot, etaInt, velBar );
        bt->Solve( NULL, G );
        VecDestroy( G );
	sw->Solve( NULL, NULL, NULL );

	WriteXDMFHeader( step + 1 );
	WriteXDMF( vFields, 3, step + 1, time, dt );
	WriteXDMF( hFields, 1, step + 1, time, dt );
	WriteXDMFFooter( step + 1 );

        bt->Setup(2);
	for( int step_i = step + 2; step_i <= step + 40000; step_i++ ) {
		time += dt;
		cout << "time step: " << step_i << "\ttime: " << time << endl;
		/* solve second order */
                bt->AssembleGSecondOrder( &top, &bot, velTop, velBot, etaInt, velTopPrev, velBotPrev, etaInt, &G );
                bt->CalcVelBar( &top, &bot, velTop, velBot, etaInt, velBar );
                bt->CalcVelBar( &top, &bot, velTopPrev, velBotPrev, etaIntPrev, velBarPrev );
                bt->Solve( velBarPrev, G );
                VecDestroy( G );
		sw->Solve( velTopPrev, velBotPrev, etaIntPrev );
		if( step_i%dumpEvery == 0 ) {
			WriteXDMFHeader( step_i );
			WriteXDMF( vFields, 3, step_i, time, dt );
			WriteXDMF( hFields, 1, step_i, time, dt );
			WriteXDMFFooter( step_i );
		}
	}

	delete vMesh;
	delete hMesh;
	delete v1bcs;
	delete v2bcs;
	delete pbcs;
	delete hbcs;
	delete velTop;
	delete velBot;
	delete etaInt;
	delete presSurf;
	delete velTopPrev;
	delete velBotPrev;
	delete etaIntPrev;
	delete sw;
        delete bt;
        delete vbbcs;
        delete velBar;
        delete velBarPrev;

	PetscFinalize();

	return EXIT_SUCCESS;
}
