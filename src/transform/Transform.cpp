#include <iostream>
#include <fstream>
#include <cmath>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include <hdf5.h>
#include <fftw3.h>

#include "Base.h"
#include "Transform.h"

using namespace std;
using std::string;

Transform::Transform( Field* _inField, int _dof, int* _nx ) {
	bool	periodic[2] 	= { false, false };

	inField 	= _inField;
	dof		= _dof;
	nx[0] 		= _nx[0];
	nx[1] 		= _nx[1];

	mesh 		= new Mesh( "transform-mesh", nx, "linear", 1, inField->mesh->min, inField->mesh->max, periodic );
	field 		= new Field( "transform-field", mesh, 1, NULL );

	grid 		= new double[nx[0]*(nx[1]+1)];
	tran		= (fftw_complex*)fftw_malloc( (nx[0]/2+1)*(nx[1]+1)*sizeof(fftw_complex) );

	real_y  	= (double*)fftw_malloc( (nx[1]+1)*sizeof(double) );
	cosine_y	= (double*)fftw_malloc( (nx[1]+1)*sizeof(double) );
	real_x		= (double*)fftw_malloc( nx[0]*sizeof(double) );
	fourier_x	= (fftw_complex*)fftw_malloc( (nx[0]/2+1)*sizeof(fftw_complex) );

	fft_x		= fftw_plan_dft_r2c_1d( nx[0], real_x, fourier_x, FFTW_ESTIMATE );
	fft_y		= fftw_plan_r2r_1d( nx[1]+1, real_y, cosine_y, FFTW_REDFT00, FFTW_ESTIMATE );
}

Transform::~Transform() {
	delete mesh;
	delete field;
	delete[] grid;
	fftw_free( tran );
	fftw_free( real_y );
	fftw_free( cosine_y );
	fftw_free( real_x );
	fftw_free( fourier_x );
	fftw_destroy_plan( fft_x );
	fftw_destroy_plan( fft_y );
}

void Transform::Trans() {
	double*		grid_t		= new double[nx[0]*(nx[1]+1)];
	double*		cosine_trans	= new double[nx[0]*(nx[1]+1)];
	double*		cosine_trans_t	= new double[nx[0]*(nx[1]+1)];

	InterpField();
	FieldToGrid();

	/* cosine transform in y dimension */
	GridTranspose( grid, grid_t, nx[0], nx[1]+1 );
	for( int i = 0; i < nx[0]; i++ ) {
		GetRow( grid_t, real_y, i, nx[1]+1 );
		fftw_execute( fft_y );
		SetRow( cosine_trans, cosine_y, i, nx[1]+1 );
	}
	/* fourier transform in x dimension */
	GridTranspose( cosine_trans, cosine_trans_t, nx[1]+1, nx[0] );
	for( int i = 0; i < nx[1]+1; i++ ) {
		GetRow( cosine_trans_t, real_x, i, nx[0] );
		fftw_execute( fft_x );
		SetRowComplex( tran, fourier_x, i, nx[0]/2+1 );
	}
	Norm();

	delete[] grid_t;
	delete[] cosine_trans;
	delete[] cosine_trans_t;
}

void Transform::Write( string filename, double minAmp ) {
	ofstream	file;

	file.open( filename.c_str() );
	for( int i = 0; i < nx[1]+1; i++ ) {
		for( int j = 0; j < nx[0]/2+1; j++ )
			if( fabs( tran[i*(nx[0]/2+1)+j][0] ) > minAmp || 
			    i == 0 || i == nx[1] || j == 0 || j == nx[0]/2 )
				file << i << "\t" << j << "\t" << tran[i*(nx[0]/2+1)+j][0] << "\t" << tran[i*(nx[0]/2+1)+j][1] << endl;
	}
	file.close();
}

void Transform::InterpField() {
	double	f[2];

	for( int node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
		inField->InterpGlobal( field->mesh->verts[node_i], f );
		field->vals[node_i][0] = f[dof];
	}
}

void Transform::FieldToGrid() {
	int	grid_i = 0, topo[2];

	for( int node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
		field->mesh->IndexToTopo( node_i, topo );
		if( topo[0] == field->mesh->nVerts[0] - 1 ) {
			continue;
		}
		grid[grid_i++] = field->vals[node_i][0];
	}
}

void Transform::GridTranspose( double* g, double* gt, int n, int m ) {
	for( int j = 0; j < m; j++ ) {
		for( int i = 0; i < n; i++ ) {
			gt[i*m+j] = g[j*n+i];
		}
	}
}

void Transform::GetRow( double* g, double* r, int row, int n ) {
	for( int i = 0; i < n; i++ ) {
		r[i] = g[row*n+i];
	}
}

void Transform::SetRow( double* g, double* r, int row, int n ) {
	for( int i = 0; i < n; i++ ) {
		g[row*n+i] = r[i];
	}
}

void Transform::SetRowComplex( fftw_complex* g, fftw_complex* r, int row, int n ) {
	for( int i = 0; i < n; i++ ) {
		g[row*n+i][0] = r[i][0];
		g[row*n+i][1] = r[i][1];
	}
}

void Transform::Norm() {
	double	fac 	= 2.0/( nx[0]*nx[1] );

	for( int i = 0; i < (nx[0]/2+1)*(nx[1]+1); i++ ) {
		tran[i][0] *= fac;
		tran[i][1] *= fac;
	}
}
