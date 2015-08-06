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

#define K 2.0
#define L 1.0
#define OMEGA sqrt( K*K + L*L )

double u0( double* x ) { return (K/OMEGA)*cos( K*x[0] + L*x[1] ); }
double v0( double* x ) { return (L/OMEGA)*cos( K*x[0] + L*x[1] ); }
double h0( double* x ) { return           cos( K*x[0] + L*x[1] ); }

void SetAnalytic( Field* velocity, Field* pressure, double time ) {
	int node_i;
	double* x;

	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		x = velocity->mesh->verts[node_i];
		velocity->vals[node_i][0] = (K/OMEGA)*cos( K*x[0] + L*x[1] - OMEGA*time );
		velocity->vals[node_i][1] = (L/OMEGA)*cos( K*x[0] + L*x[1] - OMEGA*time );
	}

	for( node_i = 0; node_i < pressure->mesh->nVertsTotal; node_i++ ) {
		x = pressure->mesh->verts[node_i];
		pressure->vals[node_i][0] = cos( K*x[0] + L*x[1] - OMEGA*time );
	}
}

void CalcDivergenceNumeric( Field* velocity, Field* divergence ) {
	double**	dv;
	int 		node_i;

	dv = new double*[2];
	dv[0] = new double[2];
	dv[1] = new double[2];

	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		velocity->InterpDerivsGlobal( velocity->mesh->verts[node_i], dv );
		divergence->vals[node_i][0] = dv[0][0] + dv[1][1];
	}

	delete[] dv[0];
	delete[] dv[1];
	delete[] dv;
}

void CalcDivergenceAnalytic( Field* divergence, double time ) {
	int node_i;
	double* x;

	for( node_i = 0; node_i < divergence->mesh->nVertsTotal; node_i++ ) {
		x = divergence->mesh->verts[node_i];
		divergence->vals[node_i][0] = -OMEGA*sin( K*x[0] + L*x[1] - OMEGA*time );
	}
}

int main( int argc, char** argv ) {
	char 			tag[] 		= "petsc";
	int			N		= 8;
	int			nx[2]		= { 12, 12 };
	double			min[2]		= { -2.0*M_PI, -2.0*M_PI };
	double			max[2]		= { +2.0*M_PI, +2.0*M_PI };
	bool			periodic[2]	= { true, true };
	Mesh*			vMesh		= new Mesh( "vMesh", nx, "legendre", N, min, max, periodic );
	Mesh*			pMesh		= new Mesh( "pMesh", nx, "legendre", N-2, min, max, periodic );
	bool			vert[2]		= { false, false };
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
	Field*			numDiv		= new Field( "numDiv", vMesh, 1, NULL );
	Field*			analDiv		= new Field( "analDiv", vMesh, 1, NULL );
	double			time		= 0.0;
	double			dt		= 0.25*(max[0] - min[0])/(nx[0]*N);
	int			timeStep;
	int			nTimeSteps	= 100;
	int 			dumpEvery	= 10;
	Field*			fields[4];
	LinearShallowWaterEqn* 	sw;

	PetscInitialize( &argc, &argv, (char*)0, tag );

	sw = new LinearShallowWaterEqn( vel, eta, 0.0, 0.0, dt );
	sw->uNullSp = sw->vNullSp = sw->pNullSp = true;

	/* initialise */
	vel->SetICFunc( 0, u0 );
	vel->SetICFunc( 1, v0 );
	eta->SetICFunc( 0, h0 );
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
	cout << "velocity error: " << FieldError( vel, vAnal, true ) << 
		"\tpressure error: " << FieldError( eta, pAnal, true ) << endl;
	CalcDivergenceNumeric( vel, numDiv );
	CalcDivergenceAnalytic( analDiv, time );
	cout << "divergence error: " << FieldError( numDiv, analDiv, true ) << endl;
	WriteXDMFHeader( 1 );
	fields[0] = vel;
	fields[1] = vAnal;
	fields[2] = numDiv;
	fields[3] = analDiv;
	WriteXDMF( fields, 4, 1, time, dt );
	fields[0] = eta;
	fields[1] = pAnal;
	fields[2] = peta;
	WriteXDMF( fields, 3, 1, time, dt );
	WriteXDMFFooter( 1 );

	for( timeStep = 2; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "time step: " << timeStep << ",\tdt: " << dt << ",\ttime: " << time << endl;
		tvel->Copy( vel );
		teta->Copy( eta );
		sw->Solve( pvel, peta, timeStep ); /* assemble the linear operator the first time this is called only */

		if( timeStep%dumpEvery == 0 ) {
			SetAnalytic( vAnal, pAnal, time );
			cout << "velocity error: " << FieldError( vel, vAnal, true ) << 
			        "\tpressure error: " << FieldError( eta, pAnal, true ) << endl;
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

	PetscFinalize();

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
	delete numDiv;
	delete analDiv;

	return 0;
}
