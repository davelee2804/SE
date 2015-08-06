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
#include "NSEqn_KIO.h"

using namespace std;
using std::string;

void SetAnalytic( Field* velocity, Field* pressure, double nu, double t ) {
	int 	node_i;
	double 	x, y;

	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		x = velocity->mesh->verts[node_i][0];
		y = velocity->mesh->verts[node_i][1];
		velocity->vals[node_i][0] = -cos(x)*sin(y)*exp(-2.0*nu*t);
		velocity->vals[node_i][1] = +sin(x)*cos(y)*exp(-2.0*nu*t);
	}
	for( node_i = 0; node_i < pressure->mesh->nVertsTotal; node_i++ ) {
		x = pressure->mesh->verts[node_i][0];
		y = pressure->mesh->verts[node_i][1];
		pressure->vals[node_i][0] = -0.25*(cos(2.0*x) + cos(2.0*y))*exp(-4.0*nu*t);
	}
}

void CalcErrors( Field* velocity, Field* pressure, Field* analVel, Field* analPres, int timeStep_i ) {
	ofstream 	file;
	char		ts[6] 		= "00000";
	char		filename[20];
	double		uErr		= FieldError( velocity, analVel, 0, true );
	double		vErr		= FieldError( velocity, analVel, 1, true );
	double		pErr		= FieldError( pressure, analPres, 0, true );

	cout << "||error||_u: " << uErr << endl;
	cout << "||error||_v: " << vErr << endl;
	cout << "||error||_p: " << pErr << endl;

	sprintf( ts, "%.5u", timeStep_i );
	sprintf( filename, "errors.%s", ts );
	file.open( filename );
	file << timeStep_i << "\t" << uErr << "\t" << vErr << "\t" << pErr << endl;
	file.close();
}

/*void UpdateBCs( Field* velocity, Field* analVel ) {
	int	nx, ny;
	int*	bottom	= velocity->bcs->GetSide( "bottom", &nx );
	int*	top	= velocity->bcs->GetSide( "top",    &nx );
	int*	left    = velocity->bcs->GetSide( "left",   &ny );
	int*	right	= velocity->bcs->GetSide( "right",  &ny );

	for( int i = 0; i < nx; i++ ) {
		velocity->vals[bottom[i]][1] = analVel
	}
}*/

int main( int argc, char** argv ) {
	char			tag[]		= "petsc";
	bool			periodic[2]	= { false, false };
	int 			nx[2]		= { 12, 12 };
	double			min[2]		= { -0.5*M_PI, -0.5*M_PI };
	double			max[2]		= { +0.5*M_PI, +0.5*M_PI };
	Mesh*			vMesh		= new Mesh( "vMesh", nx, "legendre", 7, min, max, periodic );
	Mesh*			pMesh		= new Mesh( "pMesh", nx, "legendre", 7, min, max, periodic );
	bool			vert[2]		= { false, true };
	bool			horiz[2]	= { true, false };
	//bool			vert[2]		= { true, true };
	//bool			horiz[2]	= { true, true };
	BCs*			vbcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*			pbcs		= new BCs( false, false, false, false, pMesh, 0 );
	Field*			velocity	= new Field( "velocity", vMesh, 2, vbcs );
	Field*			prevVel		= new Field( "prevVel", vMesh, 2, vbcs );
	Field*			tempVel		= new Field( "tempVel", vMesh, 2, vbcs );
	Field*			analVel		= new Field( "analVel", vMesh, 2, NULL );
	Field*			pressure	= new Field( "pressure", pMesh, 1, pbcs );
	Field*			analPres	= new Field( "analPres", pMesh, 1, NULL );
	double			nu	 	= 1.0;
	double			time		= 0.0;
	double			dx		= (max[0] - min[0])/nx[0]/8;
	double			dt1		= 0.5*(max[0] - min[0])/(2.0*nx[0]);
	double			dt2		= 0.5*dx*dx/nu;
	double			dt		= ( dt1 < dt2 ) ? dt1 : dt2;
	int			timeStep_i	= 0;
	Field*			fields[4];
	NavierStokesEqn*	ns;
	//NSEqn_KIO*		ns;

	PetscInitialize( &argc, &argv, (char)0, tag );

	cout << "dt1:\t" << dt1 << "\tdt2:\t" << dt2 << "\tdt:\t" << dt << endl;

	fields[0] = velocity;
	fields[1] = analVel;
	fields[2] = pressure;
	fields[3] = analPres;
        vMesh->Save();
	pMesh->Save();

	SetAnalytic( prevVel, pressure, nu, time );
	time += dt;
	SetAnalytic( velocity, pressure, nu, time );

	ns = new NavierStokesEqn( velocity, pressure, dt );
	//ns = new NSEqn_KIO( velocity, pressure, dt, nu );
	//ns->Setup( 2 );
	for( timeStep_i = 2; timeStep_i <= 10; timeStep_i++ ) {
		time += dt;
		SetAnalytic( analVel, analPres, nu, time );
		cout << " solving for time step " << timeStep_i << "\ttime: " << time << "\tdt: " << dt << endl;

		//UpdateBCs( velocity, analVel );

		tempVel->Copy( velocity );
		ns->Solve( velocity, pressure, prevVel, nu, 0 );
		//ns->Solve( prevVel );
		prevVel->Copy( tempVel );

		CalcErrors( velocity, pressure, analVel, analPres, timeStep_i );
		WriteXDMFHeader( timeStep_i );
		WriteXDMF( fields, 4, timeStep_i, time, dt );
		WriteXDMFFooter( timeStep_i );	
	}

	delete ns;
	delete analPres;
	delete pressure;
	delete analVel;
	delete tempVel;
	delete prevVel;
	delete velocity;
	delete pbcs;
	delete vbcs;
	delete vMesh;
	delete pMesh;

	PetscFinalize();

	return 0;
}
