#include <cmath>
#include <iostream>
#include <fstream>

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
Benjamin-Feir instability for a barotropic flow

parameters:
	L		1.0e+5  m
	H		1.0e+3  m
	U		0.1 	m/s
	f_0		5.0e-5  1/s
	beta		4.0e-11 1/ms
	g'		0.025  	m/s^2
*/

#define Ro 0.02		// Rossby number, U/f/L
#define F 1.0		// L^2/Ld^2

#define _EPS	0.02
#define _DELTA	0.02
#define _USQ	-M_PI
#define _K	1
#define _OMEGA1	-_K*1.1
#define _ALPHA	0.67876
#define _GAMMA	0.83608
#define _KAPPA	2.0
#define _B0	_KAPPA*sqrt( _ALPHA/_GAMMA )
#define _LAMBDA _GAMMA*_B0*_B0 // must be complex!
#define _OMEGA2	-_GAMMA*_B0*_B0

double* Av;
double* Ap;
double* dAv;

double LinearShear( double* x ) {
	return x[1];
}

void GenFields( Field* vel, Field* phi ) {
	ifstream	file;
	int		topo[2];

	Av  = new double[vel->mesh->nVerts[1]];
	Ap  = new double[phi->mesh->nVerts[1]];
	dAv = new double[vel->mesh->nVerts[1]];

	file.open( "Av.txt" );
	for( int y_i = 0; y_i < vel->mesh->nVerts[1]; y_i++ ) {
		file >> Av[y_i];
	}
	file.close();
	file.open( "Ap.txt" );
	for( int y_i = 0; y_i < phi->mesh->nVerts[1]; y_i++ ) {
		file >> Ap[y_i];
	}
	file.close();
	file.open( "dAv.txt" );
	for( int y_i = 0; y_i < vel->mesh->nVerts[1]; y_i++ ) {
		file >> dAv[y_i];
	}
	file.close();

	for( int vert_i = 0; vert_i < vel->mesh->nVertsTotal; vert_i++ ) {
		vel->mesh->IndexToTopo( vert_i, topo );
		vel->vals[vert_i][0] = vel->mesh->verts[vert_i][1]/M_PI - 2.0*_EPS*dAv[topo[1]]*_B0*cos( _K*vel->mesh->verts[vert_i][0] );
		vel->vals[vert_i][1] = -_K*2.0*_EPS*Av[topo[1]]*_B0*sin( _K*vel->mesh->verts[vert_i][0] );
	}
	for( int vert_i = 0; vert_i < phi->mesh->nVertsTotal; vert_i++ ) {
		phi->mesh->IndexToTopo( vert_i, topo );
		phi->vals[vert_i][0] = _USQ + 2.0*_EPS*Ap[topo[1]]*_B0*cos( _K*phi->mesh->verts[vert_i][0] );
	}
}

void CalcAnalytic( Field* vel, Field* phi, double time ) {
	int	topo[2];

	for( int vert_i = 0; vert_i < vel->mesh->nVertsTotal; vert_i++ ) {
		vel->mesh->IndexToTopo( vert_i, topo );
		vel->vals[vert_i][0] = vel->mesh->verts[vert_i][1] - 
				       2.0*_EPS*dAv[topo[1]]*_B0*cos( _K*vel->mesh->verts[vert_i][0] - _OMEGA1*time );
		vel->vals[vert_i][1] = -_K*2.0*_EPS*Av[topo[1]]*_B0*sin( _K*vel->mesh->verts[vert_i][0] - _OMEGA1*time );
	}
	for( int vert_i = 0; vert_i < phi->mesh->nVertsTotal; vert_i++ ) {
		phi->mesh->IndexToTopo( vert_i, topo );
		phi->vals[vert_i][0] = _USQ + 2.0*_EPS*Ap[topo[1]]*_B0*cos( _K*phi->mesh->verts[vert_i][0] - _OMEGA1*time );
	}
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
	int			nx[2]		= { 36, 12 };
	bool			periodic[2]	= { true, false };
	bool			noSlip[2]	= { true, true };
	bool			slip[2]		= { false, true };
	bool			nobcs[2]	= { false, false };
	double			min[2]		= { 0.0, -M_PI };
	double			max[2]		= { 6*M_PI, +M_PI };
	Mesh*			vMesh		= new Mesh( "vMesh", nx, "legendre", 8, min, max, periodic );
	Mesh*			pMesh		= new Mesh( "pMesh", nx, "legendre", 6, min, max, periodic );
	BCs*			vbcs		= new BCs( slip, slip, nobcs, nobcs, vMesh, 0 );
	BCs*			pbcs		= new BCs( false, false, false, false, pMesh, 2*vMesh->nVertsTotal - vbcs->size[0] - vbcs->size[1] );
	Field*			vel		= new Field( "vel", vMesh, 2, vbcs );
	Field*			phi		= new Field( "phi", pMesh, 1, pbcs );
	Field*			prevVel		= new Field( "prevVel", vMesh, 2, vbcs );
	Field*			prevPhi		= new Field( "prevPhi", pMesh, 1, pbcs );
	Field*			tempVel		= new Field( "tempVel", vMesh, 2, vbcs );
	Field*			tempPhi		= new Field( "tempPhi", pMesh, 1, pbcs );
	Field*			omega		= new Field( "omega", vMesh, 1, NULL );
	//double			dx		= ( vMesh->dx[0] < vMesh->dx[1] ) ? vMesh->dx[0] : vMesh->dx[1];
	//double			dt		= 0.25*dx/1.0;
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
	params.beta  = 1.9887;	// beta*L^2/U
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

	//vel->SetBCConst( "bottom", 0, +1.0 );
	//vel->SetBCConst( "top"   , 0, -1.0 );
	vel->SetBCConst( "bottom", 1,  0.0 );
	vel->SetBCConst( "top"   , 1,  0.0 );

	GenFields( vel, phi );

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
	sw->uNullSp = true;
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
	delete[] Av;
	delete[] Ap;
	delete[] dAv;

	PetscFinalize();

	return EXIT_SUCCESS;
}
