#include <iostream>
#include <fstream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include <fftw3.h>

#include "Base.h"

using namespace std;
using std::string;

void ChannelMeshCoords( Mesh* mesh ) {
	int 	index, topo[2];
	double 	x;

	for( topo[0] = 0; topo[0] < mesh->nVerts[0]; topo[0]++ ) {
		for( topo[1] = 0; topo[1] < mesh->nVerts[1]; topo[1]++ ) {
			mesh->TopoToIndex( topo, &index );
			x = M_PI - M_PI*( (double)topo[1] )/( mesh->nVerts[1] - 1 );
			mesh->verts[index][1] = mesh->min[1] + ( mesh->max[1] - mesh->min[1] )*0.5*( cos( x ) + 1.0 );
		}
	}
}

void InterpField( Field* field, Field* cField ) {
	for( int node_i = 0; node_i < cField->mesh->nVertsTotal; node_i++ ) {
		field->InterpGlobal( cField->mesh->verts[node_i], cField->vals[node_i] );
	}
}

void FieldToGrid( Field* field, double* grid ) {
	int grid_i = 0, topo[2];

	for( int node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
		field->mesh->IndexToTopo( node_i, topo );
		if( topo[0] == field->mesh->nVerts[0] - 1 ) {
			continue;
		}
		//if( topo[1] == field->mesh->nVerts[1] - 1 ) { 
		//	continue; 
		//}
		grid[grid_i++] = field->vals[node_i][0];
	}
}

void GridTranspose( double* a, double* at, int nx, int ny ) {
	for( int j = 0; j < ny; j++ ) {
		for( int i = 0; i < nx; i++ ) {
			at[i*ny+j] = a[j*nx+i];
		}
	}
}

void GetRow( double* a, double* b, int row, int n ) {
	for( int i = 0; i < n; i++ ) {
		b[i] = a[row*n+i];
	}
}

void SetRow( double* a, double* b, int row, int n ) {
	for( int i = 0; i < n; i++ ) {
		a[row*n+i] = b[i];
	}
}

void SetRowComplex( fftw_complex* a, fftw_complex* b, int row, int n ) {
	for( int i = 0; i < n; i++ ) {
		a[row*n+i][0] = b[i][0];
		a[row*n+i][1] = b[i][1];
	}
}

void WriteOutput( string filename, fftw_complex* a, int nx, int ny ) {
	ofstream	file;

	file.open( filename.c_str() );
	for( int i = 0; i < ny; i++ ) {
		for( int j = 0; j < nx; j++ ) {
			file << i << "\t" << j << "\t" << a[i*nx+j][0] << "\t" << a[i*nx+j][1] << endl;
		}
	}
	file.close();
}

void Norm( fftw_complex* a, int nx, int ny ) {
	double 	fac	= 2.0/(nx*ny);

	for( int i = 0; i < (nx/2+1)*(ny+1); i++ ) {
		a[i][0] *= fac;
		a[i][1] *= fac;
	}
}

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
	Mesh*		cMesh		= new Mesh( "c-mesh", mx, "linear", 1, min, max, periodic );
	Field*		cField		= new Field( "cField", cMesh, 1, NULL );
	double*		grid		= new double[mx[0]*(mx[1]+1)];
	int		step		= atoi( argv[1] );
	double*		real_x		= (double*)fftw_malloc( mx[0]*sizeof(double) );
	fftw_complex*	fourier_x	= (fftw_complex*)fftw_malloc( (mx[0]/2+1)*sizeof(fftw_complex) );
	double*		real_y		= (double*)fftw_malloc( (mx[1]+1)*sizeof(double) );
	double*		cosine_y	= (double*)fftw_malloc( (mx[1]+1)*sizeof(double) );
	fftw_plan	fft_x		= fftw_plan_dft_r2c_1d( mx[0], real_x, fourier_x, FFTW_ESTIMATE );
	fftw_plan	fft_y		= fftw_plan_r2r_1d( mx[1]+1, real_y, cosine_y, FFTW_REDFT00, FFTW_ESTIMATE );
	double*		grid_t		= new double[mx[0]*(mx[1]+1)];
	double*		cosine_trans	= new double[mx[0]*(mx[1]+1)];
	double*		cosine_trans_t	= new double[mx[0]*(mx[1]+1)];
	fftw_complex*	fourier_trans	= new fftw_complex[(mx[0]/2+1)*(mx[1]+1)];

	//field->Read( step );
	field->SetICFunc( 0, f );
	mesh->Save();
	WriteXDMFHeader( 0 );
	WriteXDMF( &field, 1, 0, 0.0, 0.0 );
	WriteXDMFFooter( 0 );

	//ChannelMeshCoords( cMesh );
	InterpField( field, cField );
	FieldToGrid( cField, grid );

	cMesh->Save();
	WriteXDMFHeader( 1 );
	WriteXDMF( &cField, 1, 1, 0.0, 0.0 );
	WriteXDMFFooter( 1 );

	/* cosine transform in y */
	GridTranspose( grid, grid_t, mx[0], mx[1]+1 );
	for( int i = 0; i < mx[0]; i++ ) {
		GetRow( grid_t, real_y, i, mx[1]+1 );
		fftw_execute( fft_y );
		SetRow( cosine_trans, cosine_y, i, mx[1]+1 );
	}
	/* fourier transform in x */
	GridTranspose( cosine_trans, cosine_trans_t, mx[1]+1, mx[0] );
	for( int i = 0; i < mx[1]+1; i++ ) {
		GetRow( cosine_trans_t, real_x, i, mx[0] );
		fftw_execute( fft_x );
		SetRowComplex( fourier_trans, fourier_x, i, mx[0]/2+1 );
	}
	Norm( fourier_trans, mx[0], mx[1] );

	WriteOutput( "channel.txt", fourier_trans, (mx[0]/2+1), mx[1]+1 );

	delete mesh;
	delete field;
	delete cMesh;
	delete cField;
	delete[] grid;
	fftw_free( real_x );
	fftw_free( fourier_x );
	fftw_free( real_y );
	fftw_free( cosine_y );
	fftw_destroy_plan( fft_x );
	fftw_destroy_plan( fft_y );
	delete[] grid_t;
	delete[] cosine_trans;
	delete[] cosine_trans_t;
	delete[] fourier_trans;

	return EXIT_SUCCESS;
}
