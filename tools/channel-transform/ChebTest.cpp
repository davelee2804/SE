#include <iostream>
#include <fstream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include <fftw3.h>

#include "Base.h"
#include "ChebTrans.h"

using namespace std;
using std::string;

#define NX 8

int main( int argc, char** argv ) {
	int			nx[2]		= { NX, NX };
	double			min[2] 		= { -1.0, -1.0 };
	double			max[2]		= { +1.0, +1.0 };
	bool			periodic[2]	= { false, false };
	Mesh*			mesh		= new Mesh( "mesh", nx, "linear", 1, min, max, periodic );
	ChebTrans*		ct		= new ChebTrans( NX, mesh, 0 );
	double*			f		= new double[NX+1];
	ofstream		file;

	for( int i = 0; i <= NX; i++ ) {
		//f[i] = 3.0*ct->Tn( 8, ct->x_gl_global[i] );
		f[i] = 4.0*ct->Tn( 2, ct->x_gl[i] );
	}
	ct->Transform( f );
	file.open( "chebtest.txt" );
	for( int i = 0; i <= NX; i++ ) {
		//file << i << "\t" << ct->x_gl_global[i] << "\t" << f[i] << "\t" << ct->Tj[i] << "\n";
		file << i << "\t" << ct->x_gl[i] << "\t" << f[i] << "\t" << ct->Tj[i] << "\n";
	}
	file.close();

	delete mesh;
	delete[] f;

	return EXIT_SUCCESS;
}
