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
#include "ShallowWaterEqn.h"

using namespace std;
using std::string;

double p0( double* coord ) { return 1.0 + cos( M_PI*coord[0] ); }
double v0( double* coord ) { return sin( M_PI*coord[0] ); }

void WriteToFile( Field* field, int timeStep ) {
        char            filename[40];
        ofstream        file;
        int             i, index, topo[2];

        sprintf( filename, "%s.%.5u.txt", field->name.c_str(), timeStep );

        file.open( filename );

	topo[1] = field->mesh->nVerts[1]/2;
        for( i = 0; i < field->mesh->nVerts[0]; i++ ) {
		topo[0] = i;
		field->mesh->TopoToIndex( topo, &index );
                file << i << "\t" << field->mesh->verts[index][0] << "\t" << field->vals[index][0] << endl;
        }

        file.close();
}

int main( int argc, char** argv ) {
	char 			tag[] 		= "petsc";
	int			N		= 8;
	int			nx[2]		= { 32, 1 };
	double			min[2]		= { -1.0, -0.03125 };
	double			max[2]		= { +1.0, +0.03125 };
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
	double			time		= 0.0;
	double			dt		= 0.001953125;
	int			timeStep;
	int			nTimeSteps	= 200;
	int 			dumpEvery	= 10;
	Field*			fields[2];
	double			f0		= 0.0;
	double			beta		= 0.0;
	ShallowWaterEqn* 	sw;
	double			E0, E;

	PetscInitialize( &argc, &argv, (char*)0, tag );

	sw = new ShallowWaterEqn( vel, eta, f0, beta, dt );
	sw->uNullSp = true;

	/* initialise */
	vel->SetICFunc( 0, v0 );
	eta->SetICFunc( 0, p0 );
	E0 = ShallowWaterEnergetics( vel, eta, 1.0, 0 );
	pvel->Copy( vel );
	peta->Copy( eta );

	WriteXDMFHeader( 0 );
	fields[0] = vel;
	WriteXDMF( fields, 1, 0, time, dt );
	fields[0] = eta;
	WriteXDMF( fields, 1, 0, time, dt );
	WriteXDMFFooter( 0 );
	WriteToFile( vel, 0 );
	WriteToFile( eta, 0 );

	/* 1st time step */
	time += dt;
	cout << "time step: " << 1 << ",\tdt: " << dt << ",\ttime: " << time << endl;
	sw->Solve( NULL, NULL, 1 );
	E = ShallowWaterEnergetics( vel, eta, E0, 1 );
	WriteXDMFHeader( 1 );
	fields[0] = vel;
	WriteXDMF( fields, 1, 1, time, dt );
	fields[0] = eta;
	WriteXDMF( fields, 1, 1, time, dt );
	WriteXDMFFooter( 1 );
	WriteToFile( vel, 1 );
	WriteToFile( eta, 1 );

	for( timeStep = 2; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "time step: " << timeStep << ",\tdt: " << dt << ",\ttime: " << time << endl;
		tvel->Copy( vel );
		teta->Copy( eta );
		sw->Solve( pvel, peta, timeStep ); /* assemble the linear operator the first time this is called only */

		if( timeStep%dumpEvery == 0 ) {
			E = ShallowWaterEnergetics( vel, eta, E0, timeStep );
			WriteXDMFHeader( timeStep );
			fields[0] = vel;
			WriteXDMF( fields, 1, timeStep, time, dt );
			fields[0] = eta;
			WriteXDMF( fields, 1, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteToFile( vel, timeStep );
			WriteToFile( eta, timeStep );
		}
		pvel->Copy( tvel );
		peta->Copy( teta );
	}

	PetscFinalize();

	delete sw;
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

	return 0;
}
