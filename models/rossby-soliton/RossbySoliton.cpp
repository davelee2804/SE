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
#include "OIFS.h"
#include "Advector.h"
#include "ShallowWaterEqn.h"

using namespace std;
using std::string;

#define A 0.771
#define B 0.394
#define C -0.333333333

double u0( double* x ) { return A*B*B*0.25*( 6.0*x[1]*x[1] - 9.0 )*exp( -0.5*x[1]*x[1] )/( cosh( B*x[0] )*cosh( B*x[0] ) ); }
double v0( double* x ) { return -4.0*A*B*B*B*x[1]*tanh( B*x[0] )  *exp( -0.5*x[1]*x[1] )/( cosh( B*x[0] )*cosh( B*x[0] ) ); }
double h0( double* x ) { return A*B*B*0.25*( 6.0*x[1]*x[1] + 3.0 )*exp( -0.5*x[1]*x[1] )/( cosh( B*x[0] )*cosh( B*x[0] ) ); }

double GetFieldMaxPos( Field* field );
void WriteSolitonVel( double current, double total, double* instantVel, double* accumVel, int timeStep, int dumpEvery, double dt );

void GenAnalytic( Field* vela, Field* etaa, double time ) {
	double	*x, ey2, sch;

	for( int i = 0; i < vela->mesh->nVertsTotal; i++ ) {
		x = vela->mesh->verts[i];

		ey2 = exp( -0.5*x[1]*x[1] );
		sch = 1.0/cosh( B*(x[0] - C*time - 0.395*B*B*time) );

		vela->vals[i][0] = A*B*B*0.25*( 6.0*x[1]*x[1] - 9.0 )*ey2*sch*sch;
		vela->vals[i][1] = -4.0*A*B*B*B*x[1]*tanh( B*(x[0] - C*time - 0.395*B*B*time) )*ey2*sch*sch;
	}
	for( int i = 0; i < etaa->mesh->nVertsTotal; i++ ) {
		x = etaa->mesh->verts[i];

		ey2 = exp( -0.5*x[1]*x[1] );
		sch = 1.0/cosh( B*(x[0] - C*time - 0.395*B*B*time) );

		etaa->vals[i][0] = A*B*B*0.25*( 6.0*x[1]*x[1] + 3.0 )*ey2*sch*sch;
	}
}

int main( int argc, char** argv ) {
	char 			tag[] 		= "petsc";
	int			N		= 8;
	int			nx[2]		= { 24, 6 };
	double			min[2]		= { -8.0,  0.0 };
	double			max[2]		= { +8.0, +4.0 };
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
	Field*			avel		= new Field( "avel", vMesh, 2, vbcs );
	Field*			aeta		= new Field( "aeta", pMesh, 1, hbcs );
	double			time		= 0.0;
	double			dt		= 0.0075;
	int			timeStep;
	int			nTimeSteps	= 1000;
	int 			dumpEvery	= 1;
	Field*			vFields[2];
	Field*			pFields[2];
	ShallowWaterEqn* 	sw;
	double			prevX = 0.0, currX = 0.0;
	double*			instantVel	= new double[nTimeSteps/dumpEvery];
	double*			accumVel	= new double[nTimeSteps/dumpEvery];
	SWParams		params;

	PetscInitialize( &argc, &argv, (char*)0, tag );

	params.beta = params.rho = params.L = params.H = params.U = params.g = 1.0;
	params.f0 = params.gamma = params.tau = params.nu = params.kws = 0.0;

	vFields[0] = vel;
	vFields[1] = avel;
	pFields[0] = eta;
	pFields[1] = aeta;

	sw = new ShallowWaterEqn( vel, eta, dt, 8, NULL, &params );
	sw->uNullSp = true;
	sw->pNullSp = false;

	/* initialise */
	vel->SetICFunc( 0, u0 );
	vel->SetICFunc( 1, v0 );
	eta->SetICFunc( 0, h0 );
	vel->SetBCConst( "bottom", 1, 0.0 );
	vel->SetBCConst( "top", 1, 0.0 );
	pvel->Copy( vel );
	peta->Copy( eta );

	vMesh->Save();
	pMesh->Save();

	GenAnalytic( avel, aeta, time );
	WriteXDMFHeader( 0 );
	WriteXDMF( vFields, 2, 0, time, dt );
	WriteXDMF( pFields, 2, 0, time, dt );
	WriteXDMFFooter( 0 );

	/* 1st time step */
	time += dt;
	cout << "time step: " << 1 << ",\tdt: " << dt << ",\ttime: " << time << endl;
	sw->Assemble( 1 );
	sw->Solve( NULL, NULL );
	GenAnalytic( avel, aeta, time );
	WriteXDMFHeader( 1 );
	WriteXDMF( vFields, 2, 1, time, dt );
	WriteXDMF( pFields, 2, 1, time, dt );
	WriteXDMFFooter( 1 );

	sw->Assemble( 2 );
	for( timeStep = 2; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "time step: " << timeStep << ",\tdt: " << dt << ",\ttime: " << time << endl;
		sw->Solve( pvel, peta );

		if( timeStep%dumpEvery == 0 ) {
			prevX = currX;
			currX = GetFieldMaxPos( eta );
			WriteSolitonVel( (currX - prevX)/(dumpEvery*dt), currX/time, instantVel, accumVel, timeStep, dumpEvery, dt );
			GenAnalytic( avel, aeta, time );
			WriteXDMFHeader( timeStep );
			WriteXDMF( vFields, 2, timeStep, time, dt );
			WriteXDMF( pFields, 2, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
		}
	}

	delete sw;
	delete aeta;
	delete avel;
	delete peta;
	delete pvel;
	delete vel;
	delete eta;
	delete hbcs;
	delete vbcs;
	delete pMesh;
	delete vMesh;
	delete[] instantVel;
	delete[] accumVel;

	PetscFinalize();

	return 0;
}

double GetFieldMaxPos( Field* field ) {
	int node_i;
	double max, pos = -1.0e+99;

	max = -1.0e+99;
	for( node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
		if( field->vals[node_i][0] > max ) { 
			max = field->vals[node_i][0];
			pos = field->mesh->verts[node_i][0];
		}
	}

	return pos;
}

void WriteSolitonVel( double current, double total, double* instantVel, double* accumVel, int timeStep, int dumpEvery, double dt ) {
	char		ts[6]		= "00000";
	char		filename[40];
	ofstream 	file;
	int		i;

	instantVel[timeStep/dumpEvery-1] = current;
	accumVel[timeStep/dumpEvery-1] = total;

	sprintf( ts, "%.5u", timeStep );
	sprintf( filename, "maxVel.%s.sw", ts );

	file.open( filename );
	for( i = 0; i < timeStep/dumpEvery; i++ ) {
		file << (i+1)*dumpEvery*dt << "\t" << instantVel[i] << "\t" << accumVel[i] << endl;
	}
	file.close();
}
