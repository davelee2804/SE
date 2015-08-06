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
*/

void CalcBaroclinic( Field* top, Field* bot, Field* bar, Field* topBC, Field* botBC ) {
	for( int node_i = 0; node_i < bar->mesh->nVertsTotal; node_i++ ) {
		topBC->vals[node_i][0] = top->vals[node_i][0] - bar->vals[node_i][0];
		topBC->vals[node_i][1] = top->vals[node_i][1] - bar->vals[node_i][1];
		botBC->vals[node_i][0] = bot->vals[node_i][0] - bar->vals[node_i][0];
		botBC->vals[node_i][1] = bot->vals[node_i][1] - bar->vals[node_i][1];
	}
}

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
	BCs*		hbcs		= new BCs( false, false, false, false, hMesh, 2*vMesh->nVertsTotal - vbcs->size[0] - vbcs->size[1] );
	BCs*		pbcs		= new BCs( false, false, false, false, hMesh, 0 );
	Field*		velTop		= new Field( "vel-top", vMesh, 2, vbcs );
	Field*		velBot		= new Field( "vel-bot", vMesh, 2, vbcs );
	Field*		velBar		= new Field( "vel-bar", vMesh, 2, vbcs );
	Field*		etaInt		= new Field( "eta-int", hMesh, 1, hbcs );
	Field*		presSurf	= new Field( "pres-surf", hMesh, 1, pbcs );
	Field*		velTopPrev	= new Field( "vel-top-prev", vMesh, 2, vbcs );
	Field*		velBotPrev	= new Field( "vel-bot-prev", vMesh, 2, vbcs );
	Field*		velBarPrev	= new Field( "vel-bar-prev", vMesh, 2, vbcs );
	Field*		etaIntPrev	= new Field( "eta-int-prev", hMesh, 1, hbcs );
	Field*		velTopTemp	= new Field( "vel-top-temp", vMesh, 2, vbcs );
	Field*		velBotTemp	= new Field( "vel-bot-temp", vMesh, 2, vbcs );
	Field*		velBarTemp	= new Field( "vel-bar-temp", vMesh, 2, vbcs );
	Field*		etaIntTemp	= new Field( "eta-int-temp", hMesh, 1, hbcs );
	Field*		topography	= new Field( "bot-topo", hMesh, 1, NULL );
	Field*		velTopBC	= new Field( "vel-top-bc", vMesh, 2, NULL );
	Field*		velBotBC	= new Field( "vel-bot-bc", vMesh, 2, NULL );
	double		dt		= 0.01*(max[0] - min[0])/(nx[0]*8);
	double		time		= 0.0;
	SWParams	topParams;
	SWParams	botParams;
	Field*		vFields[5];
	Field*		hFields[2];
	BaroclinicEqn*	bc;

	PetscInitialize( &argc, &argv, (char)0, tag );

	topParams.U      = botParams.U      = 1.0;
	topParams.L      = botParams.L      = 1.0e+6;
	topParams.H      = botParams.H      = 1.0e+3;
	topParams.rho_0  = botParams.rho_0  = 1.0e+3;
	topParams.f0     = botParams.f0     = 1.0e-4;
	topParams.beta   = botParams.beta   = 1.8e-11;
	topParams.nu     = botParams.nu     = 3.0e+2;
	topParams.g      = botParams.g      = 9.8;
	topParams.H_i    = 7.0e+2;
	botParams.H_i    = 3.3e+3;
	topParams.rho    = 1.000e+3;
	botParams.rho    = 1.002e+3;
	topParams.tau    = botParams.tau    = 0.1; /* 0.05?? */
	topParams.kws    = botParams.kws    = M_PI;
	botParams.gamma  = topParams.gamma  = 5.0e-8;
	botParams.gPrime = topParams.gPrime = botParams.g*( botParams.rho - topParams.rho )/botParams.rho_0;
	botParams.Htilde = topParams.Htilde = (botParams.H_i + topParams.H_i)/topParams.H;

	vFields[0] = velTop;
	vFields[1] = velBot;
	vFields[2] = velBar;
	vFields[3] = velTopBC;
	vFields[4] = velBotBC;
	hFields[0] = etaInt;
	hFields[1] = presSurf;
	vMesh->Save();
	hMesh->Save();

	bc = new BaroclinicEqn( velTop, velBot, etaInt, velBar, presSurf, topography, &topParams, &botParams, dt );
	bc->Setup( 1 );

	/* initial guess for the top layer */
	bc->tlEqn->Solve( NULL, presSurf, etaInt );
	//bc->CalcVelBar( velTop, velBot, etaInt, velBar );
	WriteXDMFHeader( 0 );
	WriteXDMF( &velTop, 1, 0, time, dt );
	WriteXDMFFooter( 0 );

	time += dt;
	bc->Solve( NULL, NULL, NULL, NULL );
	velTopPrev->Copy( velTop );
	velBotPrev->Copy( velBot );
	etaIntPrev->Copy( etaInt );
	velBarPrev->Copy( velBar );
	CalcBaroclinic( velTop, velBot, velBar, velTopBC, velBotBC );
	WriteXDMFHeader( 1 );
	WriteXDMF( vFields, 5, 1, time, dt );
	WriteXDMF( hFields, 2, 1, time, dt );
	WriteXDMFFooter( 1 );

	bc->Setup( 2 );
	for( int step_i = 2; step_i <= 100; step_i++ ) {
		time += dt;
		velTopTemp = velTop;
		velBotTemp = velBot;
		etaIntTemp = etaInt;
		velBarTemp = velBar;
		bc->Solve( velTopPrev, velBotPrev, etaIntPrev, velBarPrev );
		velTopPrev->Copy( velTopTemp );
		velBotPrev->Copy( velBotTemp );
		etaIntPrev->Copy( etaIntTemp );
		velBarPrev->Copy( velBarPrev );
		CalcBaroclinic( velTop, velBot, velBar, velTopBC, velBotBC );
		WriteXDMFHeader( step_i );
		WriteXDMF( vFields, 5, step_i, time, dt );
		WriteXDMF( hFields, 2, step_i, time, dt );
		WriteXDMFFooter( step_i );
	}

	delete vMesh;
	delete hMesh;
	delete vbcs;
	delete hbcs;
	delete pbcs;
	delete velTop;
	delete velBot;
	delete velBar;
	delete etaInt;
	delete presSurf;
	delete velTopPrev;
	delete velBotPrev;
	delete velBarPrev;
	delete etaIntPrev;
	delete velTopTemp;
	delete velBotTemp;
	delete velBarTemp;
	delete etaIntTemp;
	delete topography;
	delete velTopBC;
	delete velBotBC;
	delete bc;

	PetscFinalize();

	return EXIT_SUCCESS;
}
