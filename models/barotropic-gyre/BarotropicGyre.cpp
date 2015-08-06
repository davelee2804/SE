#include <cmath>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SWOps.h"
#include "ShallowWaterEqn.h"

using namespace std;
using std::string;

/*
Barotropic gyre model

Reference:
	Jiang, S., Jin, F. and Ghil, M. (1995) "Multiple Equilibria, Periodic, and Aperiodic Solutions
	in a Wind-Driven, Double-Gyre Shallow-Water Model" J. Phys. Oceanogr. 25, 764-786

Parameters:
	L		1.0e+6		m	( horizontal legth )
	H		5.0e+2		m	( vertical length )
	U		1.0		m/s	( velocity )
	f_0		5.0e-5		1/s	( coriolis parameter )
	beta		2.0e-11		1/ms	( coriolis gradient )
	tau_0		0.1		kg/ms^2	( wind stress, *0.8 )
	nu		3.0e+2		m^2/s	( viscosity )
	gamma		5.0e-8		1/s	( friction )
	rho		1.022e+3	kg/m^3	( density )
	g'		3.1e-2		m/s^2	( reduced gravity )
*/

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
	Field*			prevVel		= new Field( "prevVel", vMesh, 2, vbcs );
	Field*			prevPhi		= new Field( "prevPhi", pMesh, 1, pbcs );
	Field*			tempVel		= new Field( "tempVel", vMesh, 2, vbcs );
	Field*			tempPhi		= new Field( "tempPhi", pMesh, 1, pbcs );
	//double			dt		= 0.25/(nx[0]);
	double			dt		= 0.5/(3*nx[0]);
	int			timeStep	= atoi( argv[1] );
	double			time		= timeStep*dt;
	int			dumpEvery	= 50;
	int			nTimeSteps	= 80000;
	SWParams		params;
	ShallowWaterEqn* 	sw;

	PetscInitialize( &argc, &argv, (char)0, tag );

	params.L     = 1.0;
	params.H     = 1.0;
	params.U     = 1.0;
	params.rho   = 1.0;
	params.f0    = 50.0; 	// f_0*L/U
	params.beta  = 20.0;	// beta*L^2/U
	params.g     = 15.5;	// g'*H/U^2
	params.gamma = 5.0e-2;	// gamma*L/U
	params.tau   = 0.19569471624266144*0.8; // 0.8*tau*L/rho*H*U^2
	params.kws   = M_PI;
	params.nu    = 3.0e-4;	// nu/L*U

	if( timeStep > 0 ) {
		vel->Read( timeStep );
		phi->Read( timeStep );
	}

	vel->SetBCConst( "bottom", 0, 0.0 );
	vel->SetBCConst( "top"   , 0, 0.0 );
	vel->SetBCConst( "left"  , 0, 0.0 );
	vel->SetBCConst( "right" , 0, 0.0 );
	vel->SetBCConst( "bottom", 1, 0.0 );
	vel->SetBCConst( "top"   , 1, 0.0 );
	vel->SetBCConst( "left"  , 1, 0.0 );
	vel->SetBCConst( "right" , 1, 0.0 );

	vMesh->Save();
	pMesh->Save();

	if( timeStep == 0 ) {
		WriteXDMFHeader( timeStep );
		WriteXDMF( &phi, 1, timeStep, time, dt );
		WriteXDMF( &vel, 1, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
	}

	prevVel->Copy( vel );
	prevPhi->Copy( phi );

	sw = new ShallowWaterEqn( vel, phi, dt, 6, NULL, &params );
	sw->Assemble( 1 );

	time += dt;
	timeStep++;
	sw->Solve( NULL, NULL );

	WriteXDMFHeader( timeStep );
	WriteXDMF( &phi, 1, timeStep, time, dt );
	WriteXDMF( &vel, 1, timeStep, time, dt );
	WriteXDMFFooter( timeStep );

	sw->Assemble( 2 );
	for( timeStep++; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		tempVel->Copy( vel );
		tempPhi->Copy( phi );
		sw->Solve( prevVel, prevPhi );
		prevVel->Copy( tempVel );
		prevPhi->Copy( tempPhi );
		if( timeStep%dumpEvery == 0 ) {
			WriteXDMFHeader( timeStep );
			WriteXDMF( &phi, 1, timeStep, time, dt );
			WriteXDMF( &vel, 1, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
		}
	}

	delete sw;
	delete vMesh;
	delete pMesh;
	delete vbcs;
	delete pbcs;
	delete vel;
	delete phi;
	delete prevVel;
	delete prevPhi;
	delete tempVel;
	delete tempPhi;

	PetscFinalize();

	return EXIT_SUCCESS;
}
