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
Barotropic instability within a periodic channel
Unstable when beta - U'' + FU changes sign somewhere in the domain
For linear shear U = y, -1 <= y <= +1, and F = L^2/Ld^2 = 1.25, such
that flow is unstable at y = -1/1.25

parameters:
	L		1.0e+5  m
	H		1.0e+3  m
	U		0.1 	m/s
	f_0		5.0e-5  1/s
	beta		1.0e-11 1/ms	(was using 2.e-11 /m/s before)
	g'		0.02  	m/s^2
*/

#define Ro 0.02		// Rossby number, U/f/L
#define F 1.25		// L^2/Ld^2

double LinearShear( double* x ) {
	return x[1];
}

void CalcVorticity( Field* vel, Field* omega ) {
	double** gV;

	gV = new double*[2];
	gV[0] = new double[2];
	gV[1] = new double[2];

	for( int node_i = 0; node_i < vel->mesh->nVertsTotal; node_i++ ) {
		vel->InterpDerivsGlobal( vel->mesh->verts[node_i], gV );
		omega->vals[node_i][0] = gV[1][0] - gV[0][1];
	}

	delete[] gV[0];
	delete[] gV[1];
	delete[] gV;
}

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 64, 18 };
	bool			periodic[2]	= { true, false };
	bool			noSlip[2]	= { true, true };
	bool			nobcs[2]	= { false, false };
	double			min[2]		= { -M_PI, -1.0 };
	double			max[2]		= { +M_PI, +1.0 };
	Mesh*			vMesh		= new Mesh( "vMesh", nx, "legendre", 5, min, max, periodic );
	Mesh*			pMesh		= new Mesh( "pMesh", nx, "legendre", 3, min, max, periodic );
	BCs*			vbcs		= new BCs( noSlip, noSlip, nobcs, nobcs, vMesh, 0 );
	BCs*			pbcs		= new BCs( false, false, false, false, pMesh, 2*vMesh->nVertsTotal - vbcs->size[0] - vbcs->size[1] );
	Field*			vel		= new Field( "vel", vMesh, 2, vbcs );
	Field*			phi		= new Field( "phi", pMesh, 1, pbcs );
	Field*			prevVel		= new Field( "prevVel", vMesh, 2, vbcs );
	Field*			prevPhi		= new Field( "prevPhi", pMesh, 1, pbcs );
	Field*			tempVel		= new Field( "tempVel", vMesh, 2, vbcs );
	Field*			tempPhi		= new Field( "tempPhi", pMesh, 1, pbcs );
	Field*			omega		= new Field( "omega", vMesh, 1, NULL );
	double			dt		= 0.01;
	int			timeStep	= atoi( argv[1] );
	double			time		= timeStep*dt;
	int			dumpEvery	= 1;
	int			nTimeSteps	= 10000;
	SWParams		params;
	ShallowWaterEqn* 	sw;
	Field*			fields[2];

	PetscInitialize( &argc, &argv, (char)0, tag );

	params.L     = 1.0;
	params.H     = 1.0;
	params.U     = 1.0;
	params.rho   = 1.0;
	params.f0    = 50.0; 	// f_0*L/U
	params.beta  = 1.0;	// beta*L^2/U
	params.g     = 50.0;	// f_0*L/U
	params.gamma = 0.0;
	params.tau   = 0.0;
	params.kws   = 0.0;
	params.nu    = 0.0;

	if( timeStep > 0 ) {
		vel->Read( timeStep );
		phi->Read( timeStep );
	}

	fields[0] = vel;
	fields[1] = omega;

	vel->SetBCConst( "bottom", 0, 0.0 );
	vel->SetBCConst( "top"   , 0, 0.0 );
	vel->SetBCConst( "bottom", 1, 0.0 );
	vel->SetBCConst( "top"   , 1, 0.0 );

	vel->SetICFunc( 0, LinearShear );

	vMesh->Save();
	pMesh->Save();

	if( timeStep == 0 ) {
		WriteXDMFHeader( timeStep );
		WriteXDMF( fields, 2, timeStep, time, dt );
		WriteXDMF( &phi, 1, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
	}

	prevVel->Copy( vel );
	prevPhi->Copy( phi );

	sw = new ShallowWaterEqn( vel, phi, dt, 6, NULL, &params );
	sw->qgScale = 40.0; // g'*H/(f*U*L)
	sw->Assemble( 1 );

	time += dt;
	timeStep++;
	sw->Solve( NULL, NULL );

	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 2, timeStep, time, dt );
	WriteXDMF( &phi, 1, timeStep, time, dt );
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
			CalcVorticity( vel, omega );
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 2, timeStep, time, dt );
			WriteXDMF( &phi, 1, timeStep, time, dt );
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
	delete omega;

	PetscFinalize();

	return EXIT_SUCCESS;
}
