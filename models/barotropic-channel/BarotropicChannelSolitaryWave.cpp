#include <cmath>
#include <fstream>
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
1)	shift the x coordinate in the B part of the solution with the group velocity
	such that x->x - cg.t
2)	c in B is a nonlinear correction to cg and should be small w.r.t. cg (for
	small amplitude solns): check!
3)	l should be small w.r.t. k
4)	s should be O(1) & == cg!!!
5)	freeslip bcs on v (viscosity only at O(eps)) and dirichlet on p
*/

/*
propogation of a solitary wave in a single layer channel
parameters:
	L		1.0e+5  	m
	H		1.0e+3  	m
	U		0.1 		m/s
	f_0		5.0e-5  	1/s
	beta		1.9887e-11 	1/ms
	g'		0.02  		m/s^2
*/

#define EPS 0.02
#define L0 0.0 		// solitary wave number
#define OMEGA 0.01 	// solitary wave speed - should be O(1)!!!
#define ALPHA 0.35436	// nls eqn dispersive coefficient
#define GAMMA 0.81838	// nls eqn non linear coefficient
#define A0 sqrt( 2*OMEGA/GAMMA )
#define B0 sqrt( OMEGA/ALPHA )
#define C0 0.0
#define NL 12		// no. waves with wavenumber l in domain
#define K0 0.75		// base wave number
#define CP -1.14826389894846	// phase speed
#define Ro 0.02		// Rossby number, U/f/L
#define F 0.5		// L^2/Ld^2
#define CG -0.74503	// group velocity
#define LEN (NL/K0)

double* Av;
double* Ap;
double* dAv;

void LoadYStruct( Field* vel, Field* phi ) {
	ifstream	file;

	Av  = new double[vel->mesh->nVerts[1]];
	Ap  = new double[phi->mesh->nVerts[1]];
	dAv = new double[vel->mesh->nVerts[1]];

	file.open( "Av.txt" );
	for( int y_i = 0; y_i < vel->mesh->nVerts[1]; y_i++ ) { file >> Av[y_i]; }
	file.close();

	file.open( "Ap.txt" );
	for( int y_i = 0; y_i < phi->mesh->nVerts[1]; y_i++ ) { file >> Ap[y_i]; }
	file.close();

	file.open( "dAv.txt" );
	for( int y_i = 0; y_i < vel->mesh->nVerts[1]; y_i++ ) { file >> dAv[y_i]; }
	file.close();
}

void GenAnalytic( Field* phia, Field* vela, double time ) {
	int 	topo[2];
	double	X0, X1, theta, zeta;

	for( int vert_i = 0; vert_i < phia->mesh->nVertsTotal; vert_i++ ) {
		phia->mesh->IndexToTopo( vert_i, topo );
		X0 = phia->mesh->verts[vert_i][0];
		X1 = phia->mesh->verts[vert_i][0] - CG*time;
		theta = K0*X0 + L0*X1 - (CP*K0 + OMEGA)*time;
		zeta = B0*(X1 - C0*time);

		phia->vals[vert_i][0] = EPS*Ap[topo[1]]*2.0*A0*cos( theta )/cosh( zeta );
	}
	for( int vert_i = 0; vert_i < vela->mesh->nVertsTotal; vert_i++ ) {
		vela->mesh->IndexToTopo( vert_i, topo );
		X0 = vela->mesh->verts[vert_i][0];
		X1 = vela->mesh->verts[vert_i][0] - CG*time;
		theta = K0*X0 + L0*X1 - (CP*K0 + OMEGA)*time;
		zeta = B0*(X1 - C0*time);

		vela->vals[vert_i][0] = vela->mesh->verts[vert_i][1]/M_PI - EPS*dAv[topo[1]]*2.0*A0*cos( theta )/cosh( zeta );
		vela->vals[vert_i][1] = -EPS*2.0*Av[topo[1]]*A0*( 
					(K0 + L0)*sin( theta ) + B0*cos( theta )*tanh( zeta ) )/cosh( zeta );
	}
}

void RemoveShear( Field* vel, Field* velms ) {
	double	shear;
	velms->Copy( vel );
	for( int vert_i = 0; vert_i < vel->mesh->nVertsTotal; vert_i++ ) {
		shear = vel->mesh->verts[vert_i][1]/M_PI;
		velms->vals[vert_i][0] -= shear;
	}
}

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 192, 12 };
	bool			periodic[2]	= { true, false };
	bool			freeSlip[2]	= { false, true };
	bool			nobcs[2]	= { false, false };
	double			min[2]		= { -LEN*M_PI, -M_PI };
	double			max[2]		= { +LEN*M_PI, +M_PI };
	Mesh*			vMesh		= new Mesh( "vMesh", nx, "legendre", 5, min, max, periodic );
	Mesh*			pMesh		= new Mesh( "pMesh", nx, "legendre", 3, min, max, periodic );
	BCs*			vbcs		= new BCs( freeSlip, freeSlip, nobcs, nobcs, vMesh, 0 );
	BCs*			pbcs		= new BCs( true, true, false, false, pMesh, 2*vMesh->nVertsTotal - vbcs->size[0] - vbcs->size[1] );
	Field*			vel		= new Field( "vel", vMesh, 2, vbcs );
	Field*			phi		= new Field( "phi", pMesh, 1, pbcs );
	Field*			prevVel		= new Field( "prevVel", vMesh, 2, vbcs );
	Field*			prevPhi		= new Field( "prevPhi", pMesh, 1, pbcs );
	Field*			tempVel		= new Field( "tempVel", vMesh, 2, vbcs );
	Field*			tempPhi		= new Field( "tempPhi", pMesh, 1, pbcs );
	Field*			vela		= new Field( "velAnal", vMesh, 2, NULL );
	Field*			phia		= new Field( "phiAnal", pMesh, 1, NULL );
	Field*			velms		= new Field( "vel-minusShear", vMesh, 2, NULL );
	Field*			velams		= new Field( "velAnal-minusShear", vMesh, 2, NULL );
	Field*			vFields[4];
	Field*			pFields[2];
	double			dt		= 0.001;
	int			timeStep	= atoi( argv[1] );
	double			time		= timeStep*dt;
	int			dumpEvery	= 1;
	int			nTimeSteps	= 10000;
	SWParams		params;
	ShallowWaterEqn* 	sw;

	PetscInitialize( &argc, &argv, (char)0, tag );

	params.L     = 1.0;
	params.H     = 1.0;
	params.U     = 1.0;
	params.rho   = 1.0;
	params.f0    = 50.0; 	// f_0*L/U
	params.beta  = 2.0;	// beta*L^2/U
	params.g     = 50.0;	// f_0*L/U
	params.gamma = 0.0;	// gamma*L/U
	params.tau   = 0.0;     // tau*L/rho*H*U^2
	params.kws   = 0.0;
	params.nu    = 0.0;	// nu/L*U

	LoadYStruct( vel, phi );

	cout << "A0:  " << A0 << endl;
	cout << "Len: " << LEN << endl;
	cout << "c_p: " << CP << endl;
	cout << "c_g: " << CG << endl;
	cout << "C0:  " << C0 << endl;
	cout << "B0:  " << B0 << endl;

	if( timeStep > 0 ) {
		vel->Read( timeStep );
		phi->Read( timeStep );
	}
	else {
		GenAnalytic( phi, vel, time );
	}
	vFields[0] = vel;
	vFields[1] = vela;
	vFields[2] = velms;
	vFields[3] = velams;
	pFields[0] = phi;
	pFields[1] = phia;

	vel->SetBCConst( "bottom", 1, 0.0 );
	vel->SetBCConst( "top"   , 1, 0.0 );
	phi->SetBCConst( "bottom", 0, 0.0 );
	phi->SetBCConst( "top"   , 0, 0.0 );

	vMesh->Save();
	pMesh->Save();

	if( timeStep == 0 ) {
		GenAnalytic( phia, vela, 0.0 );
		WriteXDMFHeader( timeStep );
		WriteXDMF( vFields, 4, timeStep, time, dt );
		WriteXDMF( pFields, 2, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
	}

	prevVel->Copy( vel );
	prevPhi->Copy( phi );

	sw = new ShallowWaterEqn( vel, phi, dt, 6, NULL, &params );
	sw->qgScale = 40.0; // g'*H/f_0*U*L
	sw->uNullSp = true;
	sw->pNullSp = false;
	sw->presDirichlet = true;
	//sw->Assemble( 1 );

	time += dt;
	timeStep++;
	//sw->Solve( NULL, NULL );

	GenAnalytic( phia, vela, time );
	RemoveShear( vel, velms );
	RemoveShear( vela, velams );
	vel->Copy( vela );
	phi->Copy( phia );
	WriteXDMFHeader( timeStep );
	WriteXDMF( pFields, 2, timeStep, time, dt );
	WriteXDMF( vFields, 4, timeStep, time, dt );
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
			GenAnalytic( phia, vela, time );
			RemoveShear( vel, velms );
			RemoveShear( vela, velams );
			WriteXDMFHeader( timeStep );
			WriteXDMF( pFields, 2, timeStep, time, dt );
			WriteXDMF( vFields, 4, timeStep, time, dt );
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
	delete phia;
	delete velms;
	delete velams;
	delete[] Av;
	delete[] Ap;
	delete[] dAv;

	PetscFinalize();

	return EXIT_SUCCESS;
}
