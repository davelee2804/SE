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

/* Lid driven cavity flow Navier Stokes model, with Re = 10,000
 * Reference:
 *	Xiu, D & G. E. Karniadakis (2001) A Semi-Lagrangian High-Order Method for Navier-Stokes
 *	Equations, J. Comp. Phys. 172, 658-684
 */

int main( int argc, char** argv ) {
	char		tag[]		= "petsc";
	bool		periodic[2]	= { false, false };
	int 		nx[2]		= { 12, 12 };
	double		min[2]		= { 0.0, 0.0 };
	double		max[2]		= { 1.0, 1.0 };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", 8, min, max, periodic );
	bool		vert[2]		= { true, true };
	bool		horiz[2]	= { true, true };
	BCs*		vbcs		= new BCs( vert, vert, horiz, horiz, mesh, 0 );
	BCs*		pbcs		= new BCs( false, false, false, false, mesh, 0 );
	Field*		velocity	= new Field( "velocity", mesh, 2, vbcs );
	Field*		prevVel		= new Field( "prevVel", mesh, 2, vbcs );
	Field*		tempVel		= new Field( "tempVel", mesh, 2, vbcs );
	Field*		pressure	= new Field( "pressure", mesh, 1, pbcs );
	Field*		vorticity	= new Field( "vorticity", mesh, 1, NULL );
	double		nu	 	= 0.0001; 
	double		time		= 0.0;
	double		dt		= 0.5*(max[0] - min[0])/(2.0*nx[0]);
	double		U		= 1.0;
	int		timeStep_i	= 0;
	int		dumpEvery	= 10;
	Field*		fields[3];
	NavierStokesEqn*		ns;

	PetscInitialize( &argc, &argv, (char)0, tag );

	velocity->SetBCConst( "top", 0, U );
	prevVel->Copy( velocity );

	CalcVorticity( vorticity, velocity );
	WriteXDMFHeader( timeStep_i );
	fields[0] = velocity;
	fields[1] = vorticity;
	fields[2] = pressure;
	WriteXDMF( fields, 3, timeStep_i, time, dt );
	WriteXDMFFooter( timeStep_i );

	time += dt;
	timeStep_i++;
	cout << " solving for time step " << timeStep_i << "\ttime: " << time << "\tdt: " << dt << endl;

	ns = new NavierStokesEqn( velocity, pressure, dt );
	ns->Solve( velocity, pressure, NULL, true, nu, mesh->el->N/2 );

	CalcVorticity( vorticity, velocity );
	WriteXDMFHeader( timeStep_i );
	fields[0] = velocity;
	fields[1] = vorticity;
	fields[2] = pressure;
	WriteXDMF( fields, 3, timeStep_i, time, dt );
	WriteXDMFFooter( timeStep_i );

	for( timeStep_i = 2; timeStep_i <= 4000; timeStep_i++ ) {
		time += dt;
		cout << " solving for time step " << timeStep_i << "\ttime: " << time << "\tdt: " << dt << endl;

		tempVel->Copy( velocity );

		ns->Solve( velocity, pressure, prevVel, (timeStep_i==2), nu, mesh->el->N/2 );

		prevVel->Copy( tempVel );

		if( timeStep_i%dumpEvery == 0 ) {
			CalcVorticity( vorticity, velocity );
			WriteXDMFHeader( timeStep_i );
			fields[0] = velocity;
			fields[1] = vorticity;
			fields[2] = pressure;
			WriteXDMF( fields, 3, timeStep_i, time, dt );
			WriteXDMFFooter( timeStep_i );
		}
	}

	PetscFinalize();

	delete ns;
	delete vorticity;
	delete pressure;
	delete tempVel;
	delete prevVel;
	delete velocity;
	delete pbcs;
	delete vbcs;
	delete mesh;

	return 0;
}
