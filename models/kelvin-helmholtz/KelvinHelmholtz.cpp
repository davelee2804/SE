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
#include "SVV.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "Advector.h"
#include "PoissonEqn.h"
#include "HelmholtzEqn.h"
#include "NavierStokesEqn.h"

using namespace std;
using std::string;

void CalcVorticity( Field* vorticity, Field* velocity ) {
	int node_i;
	double** gV;

	gV = new double*[2];
	gV[0] = new double[2];
	gV[1] = new double[2];

	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		velocity->InterpDerivsGlobal( velocity->mesh->verts[node_i], gV );
		vorticity->vals[node_i][0] = gV[1][0] - gV[0][1];
	}

	delete[] gV[0];
	delete[] gV[1];
	delete[] gV;
}

void SetVelocity( Field* velocity ) {
	double 	a 	= 0.125;
	Mesh*	mesh	= velocity->mesh;
	double 	mid	= 0.5*(mesh->min[1] + mesh->max[1]);
	double	U	= 1.0;
	double	y;
	int	node_i;

	for( node_i = 0; node_i < mesh->nVertsTotal; node_i++ ) {
		y = mesh->verts[node_i][1];
		if( y > mid + a )      { velocity->vals[node_i][0] = +U; }
		else if( y < mid - a ) { velocity->vals[node_i][0] = -U; }
		else                   { velocity->vals[node_i][0] = U*(y - mid)/a; }
	}
}

void InterpolateInitialFields( Field* velocity, Field* pressure, int timeStep ) {
	int			nx[2]		= { 32, 8 };
	bool			periodic[2]	= { true, false };
	double			min[2]		= { 0.0, 0.0 };
	double			max[2]		= { 4.0, 1.0 };
	Mesh*			mesh		= new Mesh( "mesh", nx, "legendre", 8, min, max, periodic );
	Field*			vel0		= new Field( "velocity", mesh, 2, NULL );
	Field*			pres0		= new Field( "pressure", mesh, 1, NULL );
	int			node_i;
	double			v[2], p;

	vel0->Read( timeStep );
	pres0->Read( timeStep );

	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		vel0->InterpGlobal( velocity->mesh->verts[node_i], v );
		velocity->vals[node_i][0] = v[0];
		velocity->vals[node_i][1] = v[1];

		pres0->InterpGlobal( velocity->mesh->verts[node_i], &p );
		pressure->vals[node_i][0] = p;
	}

	delete vel0;
	delete pres0;
	delete mesh;
}

int main( int argc, char** argv ) {
	char			tag[]		= "petsc";
	bool			periodic[2]	= { true, false };
	int 			nx[2]		= { 32, 8 };
	double			min[2]		= { 0.0, 0.0 };
	double			max[2]		= { 4.0, 1.0 };
	Mesh*			mesh		= new Mesh( "vMesh", nx, "legendre", 8, min, max, periodic );
	bool			vert[2]		= { false, true };
	bool			horiz[2]	= { false, false };
	BCs*			vbcs		= new BCs( vert, vert, horiz, horiz, mesh, 0 );
	BCs*			pbcs		= new BCs( false, false, false, false, mesh, 0 );
	Field*			velocity	= new Field( "velocity", mesh, 2, vbcs );
	Field*			prevVel		= new Field( "prevVel", mesh, 2, vbcs );
	Field*			tempVel		= new Field( "tempVel", mesh, 2, vbcs );
	Field*			pressure	= new Field( "pressure", mesh, 1, pbcs );
	Field*			vorticity	= new Field( "vorticity", mesh, 1, NULL );
	double			time		= 0.0;
	double			dt		= 0.5*(max[0] - min[0])/(2.0*nx[0]);
	double			nu		= 0.0001;
	int			dumpEvery	= 10;
	Field*			fields[3];
	NavierStokesEqn*	ns;
	int			startStep	= atoi( argv[1] );
	int			timeStep_i	= ( startStep < 1 ) ? 0 : startStep;
	int			interpStep	= 350;

	PetscInitialize( &argc, &argv, (char)0, tag );
	
	fields[0] = vorticity;
	fields[1] = velocity;
	fields[2] = pressure;
	mesh->Save();

	if( startStep == interpStep ) {
		InterpolateInitialFields( velocity, pressure, startStep );
		prevVel->Copy( velocity );
	}
	else if( startStep < 1 ) {
		SetVelocity( velocity );
		prevVel->Copy( velocity );

		CalcVorticity( vorticity, velocity );
		WriteXDMFHeader( timeStep_i );
		WriteXDMF( fields, 3, timeStep_i, time, dt );
		WriteXDMFFooter( timeStep_i );
		startStep = 0;
	}
	else {
		velocity->Read( startStep );
		pressure->Read( startStep );
		prevVel->Copy( velocity );
		time = dt*startStep;
	}

	ns = new NavierStokesEqn( velocity, pressure, dt );

	time += dt;
	timeStep_i++;
	cout << " solving for time step " << timeStep_i << "\ttime: " << time << "\tdt: " << dt << endl;
	ns->Solve( velocity, pressure, NULL, nu, mesh->el->N/2 );

	CalcVorticity( vorticity, velocity );
	WriteXDMFHeader( timeStep_i );
	fields[0] = velocity;
	fields[1] = vorticity;
	fields[2] = pressure;
	WriteXDMF( fields, 3, timeStep_i, time, dt );
	WriteXDMFFooter( timeStep_i );

	for( timeStep_i++; timeStep_i <= 2000; timeStep_i++ ) {
		time += dt;
		cout << " solving for time step " << timeStep_i << "\ttime: " << time << "\tdt: " << dt << endl;
		tempVel->Copy( velocity );
		ns->Solve( velocity, pressure, prevVel, nu, mesh->el->N/2 );
		prevVel->Copy( tempVel );

		if( timeStep_i%dumpEvery == 0 ) {
			CalcVorticity( vorticity, velocity );
			WriteXDMFHeader( timeStep_i );
			WriteXDMF( fields, 3, timeStep_i, time, dt );
			WriteXDMFFooter( timeStep_i );
			WriteXDMFTemporal( timeStep_i, dumpEvery );
		}
	}

	delete ns;
	delete vorticity;
	delete pressure;
	delete tempVel;
	delete prevVel;
	delete velocity;
	delete pbcs;
	delete vbcs;
	delete mesh;

	PetscFinalize();

	return 0;
}
