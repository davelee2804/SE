#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SW.h"
#include "2LRL.h"

using namespace std;
using std::string;

/*
two layer shallow water model for the evolution of wave packets in a periodic channel

parameters:
	H 	1000 		m
	g'	0.02 		m/s^2
	L_d	25000 		m
	f_0	1.265x10^-4	s^-1
	beta	1.5x10^-11	m^-1s^-1
	U	0.1
*/

double u0( double* x ) { return 0.5/cosh( 12.0*(x[1] - 0.5) )/cosh( 12.0*(x[1] - 0.5) ); }

#define Ro 0.031621

void CalcVorticity( Field* vel, Field* vort ) {
        double** gV;

        gV = new double*[2];
        gV[0] = new double[2];
        gV[1] = new double[2];

        for( int node_i = 0; node_i < vel->mesh->nVertsTotal; node_i++ ) {
                vel->InterpDerivsGlobal( vel->mesh->verts[node_i], gV );
                vort->vals[node_i][0] = gV[1][0] - gV[0][1];
        }

        delete[] gV[0];
        delete[] gV[1];
        delete[] gV;
}

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 60, 10 };
	bool			vert[2]		= { false, true };
	bool			horiz[2]	= { false, false };
	bool			periodic[2]	= { true, false };
	double			min[2]		= { 0.0, 0.0 };
	double			max[2]		= { 6.0, 1.0 };
	Mesh*			vMesh		= new Mesh( "v-mesh", nx, "legendre", 5, min, max, periodic );
	Mesh*			hMesh		= new Mesh( "h-mesh", nx, "legendre", 3, min, max, periodic );
	BCs*			vbcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*			hbcs		= new BCs( false, false, false, false, hMesh, 2*vMesh->nVertsTotal - vbcs->size[0] - vbcs->size[1] );
	Field*			vel		= new Field( "vel", vMesh, 2, vbcs );
	Field*			eta		= new Field( "eta", hMesh, 1, hbcs );
	Field*			vort		= new Field( "vort", vMesh, 1, NULL );
	Field*			velPrev		= new Field( "vel-prev", vMesh, 2, vbcs );
	Field*			etaPrev		= new Field( "eta-prev", hMesh, 1, hbcs );
	double			dt		= 0.005;
	int			step		= atoi( argv[1] );
	double			time		= step*dt;
	int			dumpEvery	= 50;
	SWParams		params;
	Field*			fields[2];
	ShallowWaterEqn*	sw;

	PetscInitialize( &argc, &argv, (char)0, tag );

	params.L     = 1.0;
	params.U     = 1.0;
	params.f0    = 31.625;
	params.beta  = 0.09375;
	params.nu    = 1.0e-4;
	params.gamma = 0.0;
	params.g     = 2000.0;
	params.tau   = 0.0;
	params.H     = 1.0;

	vMesh->Save();
	hMesh->Save();

	fields[0] = vel;
	fields[1] = vort;

	if( step > 0 ) {
		vel->Read( step );
		eta->Read( step );
	}
	else {
		vel->SetICFunc( 0, u0 );
	}

	CalcVorticity( vel, vort );
	WriteXDMFHeader( step );
	WriteXDMF( fields, 2, step, time, dt );
	WriteXDMF( &eta, 1, step, time, dt );
	WriteXDMFFooter( step );

	/* system setup */
	sw = new ShallowWaterEqn( vel, eta, dt, 9, NULL, &params );
	sw->uNullSp = true;

	time += dt;
	step++;

	sw->Assemble( 1 );

	/* solve first order */
	velPrev->Copy( vel );
	etaPrev->Copy( eta );
	sw->Solve( NULL, NULL );

	CalcVorticity( vel, vort );
	WriteXDMFHeader( step );
	WriteXDMF( fields, 2, step, time, dt );
	WriteXDMF( &eta, 1, step, time, dt );
	WriteXDMFFooter( step );

	sw->Assemble( 2 );
	step++;
	for( ; step <= 40000; step++ ) {
		/* solve second order */
		sw->Solve( velPrev, etaPrev );
		time += dt;
		if( step%dumpEvery == 0 ) {
			CalcVorticity( vel, vort );
			WriteXDMFHeader( step );
			WriteXDMF( fields, 2, step, time, dt );
			WriteXDMF( &eta, 1, step, time, dt );
			WriteXDMFFooter( step );
			WriteXDMFTemporal( step, dumpEvery );
		}
	}

	delete vMesh;
	delete hMesh;
	delete vbcs;
	delete hbcs;
	delete vel;
	delete eta;
	delete velPrev;
	delete etaPrev;
	delete sw;

	PetscFinalize();

	return EXIT_SUCCESS;
}
