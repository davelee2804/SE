#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include <fftw3.h>

#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Transform.h"

using namespace::std;
using std::string;

double r0( double* x ) { return ( 2.0*sin(2.0*x[0]) )*( 3.0*sin(3.0*x[1]) ); }
double r1( double* x ) { return ( 2.0*sin(2.0*x[0]) + 3.0*sin(8.0*x[0]) )*( 3.0*sin(3.0*x[1]) + 4.0*sin(6.0*x[1]) ); }

void Filter( Field* field ) {
	int topo[2];

	for( int node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
		field->mesh->IndexToTopo( node_i, topo );
		if( topo[0] > 8 || topo[1] > 8 ) { field->vals[node_i][0] = 0; }
	}
}

int main( int argc, char** argv ) {
	double		min[2]		= { -M_PI, -M_PI };
	double		max[2]		= { +M_PI, +M_PI };
	bool		periodic[2]	= { false, false };
	int		nx[2]		= { 8, 8 };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", 8, min, max, periodic );
	Field*		field		= new Field( "field", mesh, 1, NULL );
	Transform*	trans		= new Transform( field );

	mesh->Save();
	trans->fMesh->Save();

	field->SetICFunc( 0, r0 );
	WriteXDMFHeader( 0 );
	WriteXDMF( &field, 1, 0, 0.0, 0.0 );
	WriteXDMFFooter( 0 );

	trans->Interp();
	WriteXDMFHeader( 1 );
	WriteXDMF( &trans->fField, 1, 1, 0.0, 0.0 );
	WriteXDMFFooter( 1 );

	field->SetICFunc( 0, r1 );
	WriteXDMFHeader( 2 );
	WriteXDMF( &field, 1, 2, 0.0, 0.0 );
	WriteXDMFFooter( 2 );

	trans->Interp();
	WriteXDMFHeader( 3 );
	WriteXDMF( &trans->fField, 1, 3, 0.0, 0.0 );
	WriteXDMFFooter( 3 );

	trans->Forward();
	WriteXDMFHeader( 4 );
	WriteXDMF( &trans->fField, 1, 4, 0.0, 0.0 );
	WriteXDMFFooter( 4 );

	Filter( trans->fField );

	trans->Backward();
	WriteXDMFHeader( 5 );
	WriteXDMF( &field, 1, 5, 0.0, 0.0 );
	WriteXDMFFooter( 5 );

	delete mesh;
	delete field;
	delete trans;

	return 1;
}
