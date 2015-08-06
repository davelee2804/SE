#include <iostream>
#include <cmath>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include <hdf5.h>

#include "Base.h"
#include "ChebTrans.h"

using namespace std;
using std::string;

ChebTrans::ChebTrans( int _n, Mesh* mesh, int dim ) {
	n 	= _n;

	/* create the gauss-lobatto points, and map to global coodinates */
	x_gl = new double[n+1];
	x_gl_global = new double[n+1];
	for( int i = 0; i <= n; i++ ) {
		x_gl[i] = cos( M_PI*( 1.0 - ((double)i)/n ) );
		x_gl_global[i] = mesh->min[dim] + ( mesh->max[dim] - mesh->min[dim] )*( 0.5*(x_gl[i] + 1.0) );
	}
	/* create the transform matrix */
	T = new double[(n+1)*(n+1)];
	for( int r = 0; r <= n; r++ ) {
		for( int c = 0; c <= n; c++ ) {
			T[r*(n+1)+c] = Tn( c, x_gl[r] );
		}
	}
	/* get the inverse transform matrix */
	Tinv = new double[(n+1)*(n+1)];
	Inv( T, Tinv, n+1 );
	Tj = new double[n+1];
}

ChebTrans::~ChebTrans() {
	delete[] x_gl;
	delete[] x_gl_global;
	delete[] T;
	delete[] Tinv;
	delete[] Tj;
}

void ChebTrans::Transform( double* in ) {
	AXEB( Tinv, in, Tj, n+1 );
}

double ChebTrans::Tn( int m, double x ) {
	return cos( m*acos( x ) );
}
