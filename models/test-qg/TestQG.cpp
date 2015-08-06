#include <cmath>
#include <fstream>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "LinSys.h"
#include "QuasiGeostrophicEqn.h"

using namespace std;
using std::string;

#define EPS 0.02
#define F1 0.5
#define BETA 1.0
#define K0 1.0
#define CP 0.39098
#define OMEGA (CP*K0)

double* A;
double* dA;
double* d2A;

void LoadYStruct( int ny ) {
	ifstream	file;

	A   = new double[ny];
	dA  = new double[ny];
	d2A = new double[ny];

	file.open( "Anew.txt" );
	for( int y_i = 0; y_i < ny; y_i++ ) { 
		file >> A[y_i];
	}
	file.close();
	file.open( "dAnew.txt" );
	for( int y_i = 0; y_i < ny; y_i++ ) { 
		file >> dA[y_i];
	}
	file.close();
	file.open( "d2Anew.txt" );
	for( int y_i = 0; y_i < ny; y_i++ ) { 
		file >> d2A[y_i];
	}
	file.close();
}

void GenAnalytic( Field* phi, Field* omega, double time ) {
	double* x;
	int topo[2];

	for( int vert_i = 0; vert_i < phi->mesh->nVertsTotal; vert_i++ ) {
		x = phi->mesh->verts[vert_i];
		phi->mesh->IndexToTopo( vert_i, topo );

		phi->vals[vert_i][0] = EPS*A[topo[1]]*cos( K0*x[0] - OMEGA*time );

		omega->vals[vert_i][0] = EPS*(d2A[topo[1]] - K0*K0*A[topo[1]])*cos( K0*x[0] - OMEGA*time );
	}
}

void GenShear( Field* phi, Field* omega, Field* vel, bool add ) {
	double* x;
	double  sfac = ( add ) ? +1.0 : -1.0;

	for( int vert_i = 0; vert_i < phi->mesh->nVertsTotal; vert_i++ ) {
		x = phi->mesh->verts[vert_i];

		phi->vals[vert_i][0] -= sfac*( x[1] + 0.20*sin( x[1] ) );

		omega->vals[vert_i][0] += sfac*( 0.2*sin( x[1] ) );

		vel->vals[vert_i][0] += sfac*( 1.0 + 0.2*cos( x[1] ) );
	}
}

void GenVelocity( Field* vel, double time ) {
	double* x;
	int topo[2];

	for( int vert_i = 0; vert_i < vel->mesh->nVertsTotal; vert_i++ ) {
		x = vel->mesh->verts[vert_i];
		vel->mesh->IndexToTopo( vert_i, topo );

		vel->vals[vert_i][0] = -EPS*dA[topo[1]]*cos( K0*x[0] - OMEGA*time );
		vel->vals[vert_i][1] = -EPS*K0*A[topo[1]]*sin( K0*x[0] - OMEGA*time );
	}
}

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 80, 20 };
	bool			periodic[2]	= { true, false };
	double			min[2]		= { -4.0*M_PI, -M_PI };
	double			max[2]		= { +4.0*M_PI, +M_PI };
	Mesh*			mesh		= new Mesh( "mesh", nx, "legendre", 5, min, max, periodic );
	BCs*			bcs		= new BCs( true, true, false, false, mesh, 0 );
	Field*			phi		= new Field( "phi-curr", mesh, 1, bcs );
	Field*			omega		= new Field( "omega-curr", mesh, 1, NULL );
	Field*			phiPrev		= new Field( "phi-prev", mesh, 1, bcs );
	Field*			omegaPrev	= new Field( "omega-prev", mesh, 1, NULL );
	Field*			phiAnal		= new Field( "phi-anal", mesh, 1, NULL );
	Field*			omegaAnal	= new Field( "omega-anal", mesh, 1, NULL );
	Field* 			velAnal 	= new Field( "vel-anal", mesh, 2, NULL );
	Field*			fields[6];
	double			dt		= 0.01;
	int			timeStep	= atoi( argv[1] );
	double			time		= timeStep*dt;
	int			dumpEvery	= 1;
	int			nTimeSteps	= 400;
	QuasiGeostrophicEqn*	qg;

	PetscInitialize( &argc, &argv, (char)0, tag );

	LoadYStruct( mesh->nVerts[1] );

	if( timeStep > 0 ) {
		phi->Read( timeStep );
		omega->Read( timeStep );
	}
	else {
		GenAnalytic( phi, omega, time );
		GenVelocity( velAnal, time );
		GenShear( phi, omega, velAnal, true );
	}

	fields[0]  = phi;
	fields[1]  = phiAnal;
	fields[2]  = omega;
	fields[3]  = omegaAnal;
	fields[5]  = velAnal;

	mesh->Save();

	//qg = new QuasiGeostrophicEqn( omega, phi, F1, BETA, dt );
	qg = new QuasiGeostrophicEqn( phi, omega, F1, BETA, 0.0, 0.0, dt );

	if( timeStep == 0 ) {
		GenAnalytic( phiAnal, omegaAnal, time );
		GenVelocity( velAnal, time );
		fields[4]  = qg->GenVel( phi );
		GenShear( phi, omega, fields[4], false );
		WriteXDMFHeader( timeStep );
		WriteXDMF( fields, 6, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
		GenShear( phi, omega, fields[4], true );
		delete fields[4];
	}

	phiPrev->Copy( phi );
	omegaPrev->Copy( omega );

	time += dt;
	timeStep++;

	GenAnalytic( phi, omega, time );
	GenVelocity( velAnal, time );
	GenShear( phi, omega, velAnal, true );
	GenAnalytic( phiAnal, omegaAnal, time );
	fields[4]  = qg->GenVel( phi );
	GenShear( phi, omega, fields[4], false );
	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 6, timeStep, time, dt );
	WriteXDMFFooter( timeStep );
	GenShear( phi, omega, fields[4], true );
	delete fields[4];

	for( timeStep++; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		//qg->Solve( phiPrev, omegaPrev );
		qg->Solve( NULL, NULL );
		if( timeStep%dumpEvery == 0 ) {
			fields[4]  = qg->GenVel( phi );
			GenVelocity( velAnal, time );
			GenAnalytic( phiAnal, omegaAnal, time );
			GenShear( phi, omega, fields[4], false );
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 6, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
			GenShear( phi, omega, fields[4], true );
			delete fields[4];
		}
	}

	delete qg;
	delete mesh;
	delete bcs;
	delete phi;
	delete omega;
	delete phiPrev;
	delete omegaPrev;
	delete phiAnal;
	delete omegaAnal;
	delete velAnal;
	delete[] A;
	delete[] dA;
	delete[] d2A;

	PetscFinalize();

	return EXIT_SUCCESS;
}
