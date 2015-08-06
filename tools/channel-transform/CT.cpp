#include <iostream>
#include <fstream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include <fftw3.h>

#include "Base.h"
#include "Transform.h"

using namespace std;
using std::string;

//double f( double* x ) { return 3.0*cos(M_PI*1.0*x[0])/* - 3.0*sin(M_PI*2.0*x[0])*/; }
double f( double* x ) { return 4.0*cos( 2.0*M_PI*x[1] )*cos( 1.0*M_PI*x[0] ); }

int main( int argc, char** argv ) {
	int		nx[2]		= { 8, 2 };
	int		mx[2]		= { 64, 16 };
	double		min[2] 		= { 0.0, 0.0 };
	double		max[2]		= { 4.0, 1.0 };
	bool		periodic[2]	= { false, false };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", 8, min, max, periodic );
	Field*		field		= new Field( "field", mesh, 1, NULL );
	//int		step		= atoi( argv[1] );
	Transform*	trans		= new Transform( field, 0, mx );

	//field->Read( step );
	field->SetICFunc( 0, f );

	mesh->Save();
	WriteXDMFHeader( 0 );
	WriteXDMF( &field, 1, 0, 0.0, 0.0 );
	WriteXDMFFooter( 0 );

	trans->Trans();
	trans->Write( "channel.txt" );

	trans->mesh->Save();
	WriteXDMFHeader( 1 );
	WriteXDMF( &trans->field, 1, 1, 0.0, 0.0 );
	WriteXDMFFooter( 1 );

	delete mesh;
	delete field;
	delete trans;

	return EXIT_SUCCESS;
}
