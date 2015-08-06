#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include <fftw3.h>

//#include <petsc.h>
//#include <petscvec.h>
//#include <petscmat.h>

#include <Base.h>

using namespace std;
using std::string;

#define NX1 240
#define NY1  20
#define NX2 960
#define NY2  80
#define LX  480.0
#define LY   40.0
#define NK  120

void InterpField( Field* field1, Field* field2 ) {
	double a;

	for( int i = 0; i < field2->mesh->nVertsTotal; i++ ) {
		field1->InterpGlobal( field2->mesh->verts[i], &a );
		field2->vals[i][0] = a;
	}
}

void SortAmps( fftw_complex* amps, int* wavenums, int n ) {
	int	kTemp, k;
	double 	amp, aTemp[2];

	for( int i = 0; i < n; i++ ) {
		amp = 0.0;
		k = 0;
		for( int j = i; j < n; j++ ) {
			if( abs( amps[j][0] ) > abs( amp ) ) {
				amp = amps[j][0];
				k = j;
			}
		}	

		kTemp = wavenums[k];
		aTemp[0] = amps[k][0];
		aTemp[1] = amps[k][1];

		wavenums[k] = wavenums[i];
		amps[k][0] = amps[i][0];
		amps[k][1] = amps[i][1];

		wavenums[i] = kTemp;
		amps[i][0] = aTemp[0];
		amps[i][1] = aTemp[1];
	}
}

void SortWaveNums( fftw_complex* amps, int* wavenums, int n ) {
	int	kTemp, k, jmin;
	double 	aTemp[2];

	for( int i = 0; i < n; i++ ) {
		jmin = 9999999;
		k = 0;
		for( int j = i; j < n; j++ ) {
			if( wavenums[j] < jmin ) {
				jmin = wavenums[j];
				k = j;
			}
		}	

		kTemp = wavenums[k];
		aTemp[0] = amps[k][0];
		aTemp[1] = amps[k][1];

		wavenums[k] = wavenums[i];
		amps[k][0] = amps[i][0];
		amps[k][1] = amps[i][1];

		wavenums[i] = kTemp;
		amps[i][0] = aTemp[0];
		amps[i][1] = aTemp[1];
	}
}

void MaxMode( Field* field, double** k, double** a, int n ) {
	int		topo[2], index, nx = field->mesh->nVerts[0]-1;
	double* 	xData		= (double*)fftw_malloc(nx*sizeof(double));
	fftw_complex* 	kData		= (fftw_complex*)fftw_malloc((nx/2+1)*sizeof(fftw_complex));
	fftw_plan	fft		= fftw_plan_dft_r2c_1d( nx, xData, kData, FFTW_MEASURE );
	double		kScale  	= 2.0*M_PI/LX;
	double		aScale  	= 2.0/NX2;
	int*		wavenums	= new int[nx/2+1];
 	
	topo[1] = field->mesh->nVerts[1]/2+2;
	for( topo[0] = 0; topo[0] < nx; topo[0]++ ) {
		field->mesh->TopoToIndex( topo, &index );
		xData[topo[0]] = field->vals[index][0];
	}
	fftw_execute( fft );

	for( int i = 0; i < nx/2+1; i++ ) {
		wavenums[i] = i;
	}

	//SortAmps( kData, wavenums, nx/2+1 );

	for( int i = 0; i < NK; i++ ) {
		k[n][i] = kScale*wavenums[i];
		a[n][i] = aScale*kData[i][0];
	}

	fftw_free( xData );
	fftw_free( kData );
	fftw_destroy_plan( fft );
	delete[] wavenums;
}

void GetMode( Field* field, int kn, double* a, int n ) {
	int		topo[2], index, nx = field->mesh->nVerts[0]-1;
	double* 	xData		= (double*)fftw_malloc(nx*sizeof(double));
	fftw_complex* 	kData		= (fftw_complex*)fftw_malloc((nx/2+1)*sizeof(fftw_complex));
	fftw_plan	fft		= fftw_plan_dft_r2c_1d( nx, xData, kData, FFTW_MEASURE );
	double		aScale  	= 2.0/NX2;
 	
	topo[1] = field->mesh->nVerts[1]/2+2;
	for( topo[0] = 0; topo[0] < nx; topo[0]++ ) {
		field->mesh->TopoToIndex( topo, &index );
		xData[topo[0]] = field->vals[index][0];
	}
	fftw_execute( fft );

	a[n] = aScale*kData[kn][0];

	fftw_free( xData );
	fftw_free( kData );
	fftw_destroy_plan( fft );
}

void GenMaxAmp( Field* field, double* maxAmp, double* maxPos, int s ) {
	maxAmp[s] = 0.0;

	for( int i = 0; i < field->mesh->nVertsTotal; i++ ) {
		if( field->vals[i][0] > maxAmp[s] ) {
			maxAmp[s] = field->vals[i][0];
			maxPos[s] = field->mesh->verts[i][0];
		}
	}
}

int main( int argc, char** argv ) {
	int		nx1[2]		= { NX1, NY1 };
	int		nx2[2]		= { NX2, NY2 };
	double		min[2]		= { 0.0, -0.5*LY };
	double		max[2]		= { LX,  +0.5*LY };
	bool		periodic[2]	= { true, false };
	Mesh*		mesh1		= new Mesh( "mesh", nx1, "legendre", 4, min, max, periodic );
	Mesh*		mesh2		= new Mesh( "mesh-2", nx2, "linear", 1, min, max, periodic );
	int		minStep		= atoi( argv[1] );
	int		maxStep		= atoi( argv[2] );
	int 		step		= atoi( argv[3] );
	int		nSteps		= (maxStep-minStep)/step + 2;
	string		name		= argv[4];
	int		kn		= atoi( argv[5] );
	Field*		field1		= new Field( name.c_str(), mesh1, 1, NULL );
	Field*		field2		= new Field( "field-2",    mesh2, 1, NULL );
	int 		i, j;
	double**	k		= new double*[nSteps];
	double**	a		= new double*[nSteps];
	double*		amp		= new double[nSteps];
	double*		maxAmp		= new double[nSteps];
	double*		maxPos		= new double[nSteps];
	ofstream	file;
	char		filename[50];

	for( i = 0; i < nSteps; i++ ) {
		k[i] = new double[NK];
		a[i] = new double[NK];
	}

	cout << "generating growth rates for " << name.c_str() << endl;

	i = 0;
	for( int s = minStep; s <= maxStep; s += step ) {
		cout << "step: " << s << endl;
		field1->Read( s );
		InterpField( field1, field2 );
		MaxMode( field2, k, a, i );
		GetMode( field2, kn, amp, i );
		cout << "k max: " << LX/2.0/M_PI*k[i][0] << "\tamp: " << a[i][0] << "\tk max: " << LX/2.0/M_PI*k[i][1] << "\tamp: " << a[i][1] << endl;
		GenMaxAmp( field1, maxAmp, maxPos, i );
		i++;
	}
	mesh2->Save();

	file.open( "k.dat" );
	for( i = 0; i < nSteps; i++ ) {
		for( j = 0; j < NK; j++ ) {
			file << k[i][j] << "\t";
		}
		file << endl;
	}
	file.close();

	file.open( "a.dat" );
	for( i = 0; i < nSteps; i++ ) {
		for( j = 0; j < NK; j++ ) {
			file << a[i][j] << "\t";
		}
		file << endl;
	}
	file.close();

	file.open( "cg.dat" );
	for( i = 0; i < nSteps; i++ ) {
		if( i == 0 ) {
			file << maxPos[i] << "\t" << maxAmp[i] << "\t0.0" << endl;
		}
		else {
			file << maxPos[i] << "\t" << maxAmp[i] << "\t" << maxPos[i]-maxPos[i-1] << endl;
		}
	}
	file.close();

	sprintf( filename, "amp_%.2u.dat", kn );
	file.open( filename );
	for( i = 0; i < nSteps; i++ ) {
		file << amp[i] << endl;
	}
	file.close();

	delete field1;
	delete field2;
	delete mesh1;
	delete mesh2;
	for( i = 0; i < nSteps; i++ ) {
		delete[] k[i];
		delete[] a[i];
	}
	delete[] k;
	delete[] a;
	delete[] amp;
	delete[] maxAmp;
	delete[] maxPos;

	return EXIT_SUCCESS;
}
