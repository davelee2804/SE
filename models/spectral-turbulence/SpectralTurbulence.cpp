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
#include "NavierStokesSpectral.h"

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

void InitVelocity( Field* velocity, NavierStokesSpectral* ns ) {
	int		nx		= velocity->mesh->nVerts[0] - 1;
	int		ny		= velocity->mesh->nVerts[1] - 1;
	fftw_complex*	fourier		= (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	double*		real		= (double*)fftw_malloc( nx*ny*sizeof(double) );
	int		mode_i;
	int		kx, ky, l;
	fftw_plan	backward;

	for( mode_i = 0; mode_i < nx*(ny/2+1); mode_i++ ) {
		fourier[mode_i][0] = 0.0;
		fourier[mode_i][1] = 0.0;
	}
	srand( 101 );
	for( l = 0; l < 4; l++ ) {
		do {
			kx = rand()%8;
			ky = rand()%8 - 4;
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
	ns->MapFromArray( velocity, real, 0 );
	//ns->MapFromArray( velocity, real, 1 );

	fftw_destroy_plan( backward );
	fftw_free( fourier );
	fftw_free( real );
}

int main( int argc, char** argv ) {
	int			nx[2]		= { 256, 256 };
	bool			periodic[2]	= { true, true };
	double 			min[2]		= { 0.0, 0.0 };
	double			max[2]		= { 1.0, 1.0 };
	Mesh*			mesh		= new Mesh( "mesh", nx, "quadratic", 2, min, max, periodic );
	Field*			velocity	= new Field( "velocity", mesh, 2, NULL );
	Field*			pressure	= new Field( "pressure", mesh, 1, NULL );
	Field*  		vorticity	= new Field( "vorticity", mesh, 1, NULL );
	Field*			velPrev		= new Field( "velPrev", mesh, 2, NULL );
	Field*			tempVel		= new Field( "tempVel", mesh, 2, NULL );
	Field*			fields[3];
	int 			timeStep	= 0;
	int			dumpEvery	= 10;
	double			time		= 0.0;
	double			dt		= 0.5*1.0/256.0;
	NavierStokesSpectral*	ns		= new NavierStokesSpectral( velocity, pressure, 0.0 );

	fields[0] = velocity;
	fields[1] = pressure;
	fields[2] = vorticity;

	InitVelocity( velocity, ns );
	velPrev->Copy( velocity );
	CalcVorticity( vorticity, velocity );
	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 3, 0, 0.0, dt );
	WriteXDMFFooter( 0 );

	time += dt;
	timeStep++;
	cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
	ns->Solve( dt, NULL );
	if( timeStep%dumpEvery == 0 ) {
		CalcVorticity( vorticity, velocity );
		WriteXDMFHeader( timeStep );
		WriteXDMF( fields, 3, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
	}

	for( timeStep = 2; timeStep <= 10000; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		tempVel->Copy( velocity );
		ns->Solve( dt, velPrev );
		velPrev->Copy( tempVel );
		if( timeStep%dumpEvery == 0 ) {
			CalcVorticity( vorticity, velocity );
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 3, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
		}
	}

	delete ns;
	delete velocity;
	delete pressure;
	delete vorticity;
	delete velPrev;
	delete tempVel;
	delete mesh;

	return 0;
}
