#include <iostream>
#include <string>
#include <fstream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "QuadPoint.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "SVV.h"
#include "LinearShallowWaterEqn.h"

using namespace std;
using std::string;

#define C 1.0
#define F0 0.0
#define BETA 1.0
#define K 2.0
#define OMEGA 0.5*(C*K + sqrt(C*C*K*K + 4.0*C*BETA))

double u0( double* x ) { return ( ( BETA*x[1]*exp( -0.5*BETA*x[1]*x[1]/C ) )/( OMEGA - C*K ) )*cos( K*x[0] ); }
double v0( double* x ) { return exp( -0.5*BETA*x[1]*x[1]/C )*sin( K*x[0] ); }
double h0( double* x ) { return ( ( BETA*x[1]*exp( -0.5*BETA*x[1]*x[1]/C ) )/( C*OMEGA - C*C*K ) )*cos( K*x[0] ); }

void SetAnalytic( Field* velocity, Field* pressure, double time ) {
	int node_i;
	double* x;

	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		x = velocity->mesh->verts[node_i];
		velocity->vals[node_i][0] = ( ( BETA*x[1]*exp( -0.5*BETA*x[1]*x[1]/C ) )/( OMEGA - C*K ) )*cos( K*x[0] - OMEGA*time );
		velocity->vals[node_i][1] = exp( -0.5*BETA*x[1]*x[1]/C )*sin( K*x[0] - OMEGA*time );
	}

	for( node_i = 0; node_i < pressure->mesh->nVertsTotal; node_i++ ) {
		x = pressure->mesh->verts[node_i];
		pressure->vals[node_i][0] = ( ( BETA*x[1]*exp( -0.5*BETA*x[1]*x[1]/C ) )/( C*OMEGA - C*C*K ) )*cos( K*x[0] - OMEGA*time );
	}
}

int main( int argc, char** argv ) {
	char 			tag[] 		= "petsc";
	int			N		= 8;
	int			nx[2]		= { 24, 12 };
	double			min[2]		= { -2.0*M_PI, -M_PI };
	double			max[2]		= { +2.0*M_PI, +M_PI };
	bool			periodic[2]	= { true, false };
	Mesh*			vMesh		= new Mesh( "vMesh", nx, "legendre", N, min, max, periodic );
	Mesh*			pMesh		= new Mesh( "pMesh", nx, "legendre", N-2, min, max, periodic );
	bool			vert[2]		= { false, true };
	bool			horiz[2]	= { false, false };
	BCs*			vbcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*			hbcs		= new BCs( false, false, false, false, pMesh, 2*vMesh->nVertsTotal - vbcs->size[0] - vbcs->size[1] );
	Field*			vel		= new Field( "vel", vMesh, 2, vbcs );
	Field*			eta		= new Field( "eta", pMesh, 1, hbcs );
	Field*			pvel		= new Field( "pvel", vMesh, 2, vbcs );
	Field*			peta		= new Field( "peta", pMesh, 1, hbcs );
	Field*			tvel		= new Field( "tvel", vMesh, 2, vbcs );
	Field*			teta		= new Field( "teta", pMesh, 1, hbcs );
	Field*			vAnal		= new Field( "vAnal", vMesh, 2, NULL );
	Field*			pAnal		= new Field( "pAnal", pMesh, 1, NULL );
	double			time		= 0.0;
	double			dt		= 0.25*(max[0] - min[0])/(nx[0]*N);
	int			timeStep;
	int			nTimeSteps	= 100;
	double			vxErr[100];
	double			vyErr[100];
	double			pErr[100];
	int 			dumpEvery	= 1;
	Field*			fields[2];
	LinearShallowWaterEqn* 	sw;
	ofstream		file;

	PetscInitialize( &argc, &argv, (char*)0, tag );

	sw = new LinearShallowWaterEqn( vel, eta, F0, BETA, 1.0, 0.0, 0.0, 0.0, dt, 0, NULL );
	sw->uNullSp = true;

	vMesh->Save();
	pMesh->Save();

	vxErr[0] = vyErr[0] = pErr[0] = 0.0;

	/* initialise */
	vel->SetICFunc( 0, u0 );
	vel->SetICFunc( 1, v0 );
	eta->SetICFunc( 0, h0 );
	vel->SetBCConst( "bottom", 1, 0.0 );
	vel->SetBCConst( "top", 1, 0.0 );
	pvel->Copy( vel );
	peta->Copy( eta );

	SetAnalytic( vAnal, pAnal, time );
	WriteXDMFHeader( 0 );
	fields[0] = vel;
	fields[1] = vAnal;
	WriteXDMF( fields, 2, 0, time, dt );
	fields[0] = eta;
	fields[1] = pAnal;
	WriteXDMF( fields, 2, 0, time, dt );
	WriteXDMFFooter( 0 );

	/* 1st time step */
	time += dt;
	cout << "time step: " << 1 << ",\tdt: " << dt << ",\ttime: " << time << endl;
	sw->Solve( NULL, NULL, 1 );

	SetAnalytic( vAnal, pAnal, time );
	vxErr[timeStep] = FieldError( vel, vAnal, 0, true );
	vyErr[timeStep] = FieldError( vel, vAnal, 1, true );
	pErr[timeStep]  = FieldError( eta, pAnal, 0, true );
	cout << "vx error: " << vxErr[timeStep] << "\tvx error: " << vyErr[timeStep] << "\tp  error: " << pErr[timeStep] << endl;

	WriteXDMFHeader( 1 );
	fields[0] = vel;
	fields[1] = vAnal;
	WriteXDMF( fields, 2, 1, time, dt );
	fields[0] = eta;
	fields[1] = pAnal;
	WriteXDMF( fields, 2, 1, time, dt );
	WriteXDMFFooter( 1 );

	for( timeStep = 2; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "time step: " << timeStep << ",\tdt: " << dt << ",\ttime: " << time << endl;
		tvel->Copy( vel );
		teta->Copy( eta );
		sw->Solve( pvel, peta, timeStep ); /* assemble the linear operator the first time this is called only */

		if( timeStep%dumpEvery == 0 ) {
			SetAnalytic( vAnal, pAnal, time );
			vxErr[timeStep] = FieldError( vel, vAnal, 0, true );
			vyErr[timeStep] = FieldError( vel, vAnal, 1, true );
			pErr[timeStep]  = FieldError( eta, pAnal, 0, true );
			cout << "vx error: " << vxErr[timeStep] << "\tvx error: " << vyErr[timeStep] << "\tp  error: " << pErr[timeStep] << endl;
			WriteXDMFHeader( timeStep );
			fields[0] = vel;
			fields[1] = vAnal;
			WriteXDMF( fields, 2, timeStep, time, dt );
			fields[0] = eta;
			fields[1] = pAnal;
			WriteXDMF( fields, 2, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
		}
		pvel->Copy( tvel );
		peta->Copy( teta );
	}

	file.open( "rk_err.txt" );
	for( int i = 0; i < nTimeSteps; i++ ) {
		file << i << "\t" << i*dt << "\t" << vxErr[i] << "\t" << vyErr[i] << "\t" << pErr[i] << endl;
	}
	file.close();

	delete sw;
	delete vAnal;
	delete pAnal;
	delete teta;
	delete tvel;
	delete peta;
	delete pvel;
	delete vel;
	delete eta;
	delete hbcs;
	delete vbcs;
	delete pMesh;
	delete vMesh;

	PetscFinalize();

	return 0;
}
