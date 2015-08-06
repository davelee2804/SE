#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

#include <hdf5.h>
#include <fftw3.h>

#include "QuadPoint.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Advector.h"
#include "QuasiGeostrophicSpectral.h"

using namespace std;
using std::string;

void InitFields( Field* vorticity, Field* velocity, QuasiGeostrophicSpectral* gq ) {
	int		nx		= velocity->mesh->nVerts[0] - 1;
	int		ny		= velocity->mesh->nVerts[1] - 1;
	fftw_complex*	fourier		= (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	double*		real		= (double*)fftw_malloc( nx*ny*sizeof(double) );
	int		kx, ky, l;
	fftw_plan	backward;
        double** 	gV;

	for( int mode_i = 0; mode_i < nx*(ny/2+1); mode_i++ ) {
		fourier[mode_i][0] = 0.0;
		fourier[mode_i][1] = 0.0;
	}
	srand( 101 );
	for( l = 0; l < 4; l++ ) {
		do {
			kx = rand()%16;
			ky = rand()%16 - 8;
			if( ky < 0 ) {
				ky = ny - ky;
			}
		} while( kx == 0 || ky == 0 );
		fourier[ky*(nx/2+1)+kx][0] = 1.0;
		//fourier[(nx-k)*(nx/2+1)+k][0] = 1.0;
		cout << "kx[" << l << "]: " << kx << "\tky[" << l << "]: " << ky << endl;
	}

	backward = fftw_plan_dft_c2r_2d( nx, ny, fourier, real, FFTW_ESTIMATE );
	fftw_execute( backward );
	gq->MapFromArray( velocity, real, 0 );

	fftw_destroy_plan( backward );
	fftw_free( fourier );
	fftw_free( real );

        gV = new double*[2];
        gV[0] = new double[2];
        gV[1] = new double[2];

        for( int node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
                velocity->InterpDerivsGlobal( velocity->mesh->verts[node_i], gV );
                vorticity->vals[node_i][0] = gV[1][0] - gV[0][1];
        }

        delete[] gV[0];
        delete[] gV[1];
        delete[] gV;
}

int main( int argc, char** argv ) {
	int				nx[2]		= { 128,  128 };
	bool				periodic[2]	= { true, true };
	double 				min[2]		= { 0.0, 0.0 };
	double				max[2]		= { 1.0, 1.0 };
	Mesh*				mesh		= new Mesh( "mesh", nx, "quadratic", 2, min, max, periodic );
	Field*				velocity	= new Field( "velocity", mesh, 2, NULL );
	Field*				streamFunc	= new Field( "streamFunc", mesh, 1, NULL );
	Field*  			vorticity	= new Field( "vorticity", mesh, 1, NULL );
	Field*				prevVort	= new Field( "prevVort", mesh, 1, NULL );
	Field*				tempVort	= new Field( "tempVort", mesh, 1, NULL );
	Field*				prevVel 	= new Field( "prevVel", mesh, 2, NULL );
	Field*				tempVel 	= new Field( "tempVel", mesh, 2, NULL );
	Field*				fields[3];
	int 				timeStep	= 0;
	int				dumpEvery	= 2;
	double				time		= 0.0;
	double				dt		= 0.5*1.0/128.0;
	QuasiGeostrophicSpectral*	qg		= new QuasiGeostrophicSpectral( vorticity, streamFunc, velocity, 200.0 );

	fields[0] = vorticity;
	fields[1] = streamFunc;
	fields[2] = velocity;

	InitFields( vorticity, velocity, qg );
	mesh->Save();
	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 1, 0, 0.0, dt );
	WriteXDMFFooter( 0 );

	/* first time step */
	time += dt;
	timeStep++;
	cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
	qg->Solve( dt, NULL, NULL );
	prevVort->Copy( vorticity );
	prevVel->Copy( velocity );
	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 3, timeStep, time, dt );
	WriteXDMFFooter( timeStep );

	/* second time step */
	time += dt;
	timeStep++;
	cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
	qg->Solve( dt, NULL, NULL );
	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 3, timeStep, time, dt );
	WriteXDMFFooter( timeStep );

	for( timeStep = 3; timeStep <= 1000; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		tempVort->Copy( vorticity );
		tempVel->Copy( velocity );

		qg->Solve( dt, prevVort, prevVel );

		prevVort->Copy( tempVort );
		prevVel->Copy( tempVel );
		if( timeStep%dumpEvery == 0 ) {
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 3, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
		}
	}

	delete qg;
	delete velocity;
	delete streamFunc;
	delete vorticity;
	delete prevVort;
	delete tempVort;
	delete prevVel;
	delete tempVel;
	delete mesh;

	return 0;
}
