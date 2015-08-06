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

double CalcTotalDivergence( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* etaInt ) {
	double	**du1, **du2, h1, h2, e, div = 0.0, div_e = 0.0, *coord, detJac, weight, gCoord[2];

	du1 = new double*[2];
	du2 = new double*[2];
	du1[0] = new double[2];
	du1[1] = new double[2];
	du2[0] = new double[2];
	du2[1] = new double[2];

	for( int el_i = 0; el_i < velTop->mesh->nElsTotal; el_i++ ) {
		div_e = 0.0;
		for( int pt_i = 0; pt_i < velTop->mesh->el->nPoints; pt_i++ ) {
			coord  = velTop->mesh->el->quadPts[pt_i]->coord;
			weight = velTop->mesh->el->quadPts[pt_i]->weight;
			detJac = velTop->mesh->DetJac( el_i, pt_i );
			etaInt->InterpLocal( el_i, coord, &e );
			velTop->mesh->LocalToGlobal( coord, el_i, gCoord );
			velTop->InterpDerivsGlobal( gCoord, du1 );
			velBot->InterpDerivsGlobal( gCoord, du2 );
			h1 = tp->H - e;
			h2 = bp->H + e;
			div_e += detJac*weight*( h1*( du1[0][0] + du1[1][1] ) + h2*( du2[0][0] + du2[1][1] ) )/( tp->H + bp->H );
		}
		div += div_e;
	}

	delete[] du1[0];
	delete[] du1[1];
	delete[] du2[0];
	delete[] du2[1];
	delete[] du1;
	delete[] du2;

	return div_e;
}

int main( int argc, char** argv ) {
	char		tag[6]		= "petsc";
	int		nx[2]		= { 32, 64 };
	bool		noSlip[2]	= { true, true };
	bool		periodic[2]	= { false, false };
	double		min[2]		= { 0.0, 0.0 };
	double		max[2]		= { 1.0, 2.0 };
	Mesh*		vMesh		= new Mesh( "v-mesh", nx, "quadratic", 2, min, max, periodic );
	Mesh*		hMesh		= new Mesh( "h-mesh", nx, "linear", 1, min, max, periodic );
	BCs*		vbbcs		= new BCs( noSlip, noSlip, noSlip, noSlip, vMesh, 0 );
	BCs*		v1bcs		= new BCs( noSlip, noSlip, noSlip, noSlip, vMesh, 0 );
	BCs*		v2bcs		= new BCs( noSlip, noSlip, noSlip, noSlip, vMesh, 2*vMesh->nVertsTotal - v1bcs->size[0] - v1bcs->size[1] );
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
	double		dt		= 0.25*(max[0] - min[0])/(nx[0]*4);
	double		time		= 0.0;
	int		step		= atoi( argv[1] );
	int		dumpEvery	= 40;
	Params		topParams;
	Params		botParams;
	Params		barParams;
	Field*		vFields[4];
	Field*		hFields[1];
	MomentumEqn*	tl		= NULL;
	BarotropicEqn*	bt;
	SW2L*		sw;
	Vec		G;

	PetscInitialize( &argc, &argv, (char)0, tag );

	topParams.f0    = botParams.f0    = barParams.f0    = 1.0e+2;
	topParams.beta  = botParams.beta  = barParams.beta  = 1.8e+1;
	topParams.nu    = botParams.nu    = barParams.nu    = 3.0e-4;
	topParams.gamma = 0.0;
	botParams.gamma                   = barParams.gamma = 5.0e-2;
	topParams.pFac  = botParams.pFac  = barParams.pFac  = 1.0e+2;
	topParams.hFac  = 0.0;
	botParams.hFac                    = barParams.hFac  = 1.0e+2;
	topParams.bFac  = 0.0;
	botParams.bFac                    = barParams.bFac  = 0.0;
	topParams.tau                     = barParams.tau   = 0.1;	/* 0.05?? */
	botParams.tau   = 0.0;
	topParams.k                       = barParams.k     = M_PI;
	botParams.k     = 0.0;
	topParams.H	= 0.7;
	botParams.H 	= 3.3;
	barParams.H     = 4.0;
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
	else { /* initial guess for the top layer */
		tl = new MomentumEqn( velTop, etaInt, velTopPrev, etaIntPrev, NULL, NULL, &topParams, dt );
		tl->Assemble( 1 );
		tl->Solve( 1 );
		tl->Update();
		delete tl;
		WriteXDMFHeader( 0 );
		WriteXDMF( &velTop, 1, 0, time, dt );
		WriteXDMFFooter( 0 );
		step = 0;
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

	cout << "divergence: " << CalcTotalDivergence( &topParams, &botParams, velTop, velBot, etaInt ) << endl;

	WriteXDMFHeader( step + 1 );
	WriteXDMF( vFields, 4, step + 1, time, dt );
	WriteXDMF( hFields, 1, step + 1, time, dt );
	WriteXDMFFooter( step + 1 );

	bt->Setup( 2 );
	sw->Assemble( 2 );
	for( int step_i = step + 2; step_i <= step + 10000; step_i++ ) {
		/* solve second order */
		//bt->AssembleGSecondOrder( &topParams, &botParams, velTop, velBot, etaInt, velTopPrev, velBotPrev, etaIntPrev, &G );
	//bt->AssembleGFirstOrder( &topParams, &botParams, velTop, velBot, etaInt, &G );
	bt->AssembleGSecondOrder( &topParams, &botParams, velTop, velBot, etaInt, velTopPrev, velBotPrev, etaInt, &G );
		bt->CalcVelBar( &topParams, &botParams, velTop, velBot, etaInt, velBar );
		bt->CalcVelBar( &topParams, &botParams, velTopPrev, velBotPrev, etaIntPrev, velBarPrev );
		bt->Solve( velBarPrev, G );
		sw->Solve( 2 );
		VecDestroy( G );
		cout << "divergence: " << CalcTotalDivergence( &topParams, &botParams, velTop, velBot, etaInt ) << endl;
		time += dt;
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
	delete tl;
	delete bt;
	delete sw;

	PetscFinalize();

	return EXIT_SUCCESS;
}
