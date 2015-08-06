#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <hdf5.h>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Legendre.h"
#include "Chebyshev.h"
#include "Mesh.h"

using namespace std;
using std::string;

Mesh::Mesh( string _name, int* _nEls, string elType, int _N, double* _min, double* _max, bool* _periodic ) {
	static int dim_i, el_i, el_j, vert_i, topo[2];
	static double lCoord[2];

	name = _name;
	N = _N;
	dim = 2;

	if( elType == "chebyshev" ) {
		el = new Chebyshev( N );
	}
	else if( elType == "legendre" ) {
		el = new Legendre( N );
	}
	else if( elType == "lagrange" ) {
		el = new Lagrange( N );
	}
	else if( elType == "trig" ) {
		el = new Trig( N );
	}
	else if( elType == "linear" ) {
		N = 1;
		el = new Linear( N );
	}
	else if( elType == "quadratic" ) {
		N = 2;
		el = new Quadratic( N );
	}
	else {
		cerr << "unrecognised element type: " << elType << endl;
		exit( 1 );
	}

	for( dim_i = 0; dim_i < 2; dim_i++ ) {
		nEls[dim_i] = _nEls[dim_i];
		nVerts[dim_i] = nEls[dim_i]*N + 1;
		min[dim_i] = _min[dim_i];
		max[dim_i] = _max[dim_i];
		periodic[dim_i] = _periodic[dim_i];
		dx[dim_i] = (max[dim_i] - min[dim_i])/nEls[dim_i];
	}

	/* assign the global coordinates */
	nElsTotal = nEls[0]*nEls[1];
	nVertsTotal = nVerts[0]*nVerts[1];
	verts = new double*[nVertsTotal];
	for( vert_i = 0; vert_i < nVertsTotal; vert_i++ ) {
		verts[vert_i] = new double[2];
		IndexToTopo( vert_i, topo );
		el_i = topo[0]/N;
		el_j = topo[1]/N;
		lCoord[0] = el->xi[topo[0]%N][0];
		lCoord[1] = el->xi[topo[1]%N][0]; /* element geometry in the same in each dimension */
		verts[vert_i][0] = min[0] + el_i*dx[0] + ((lCoord[0] + 1.0)/2.0)*dx[0];
		verts[vert_i][1] = min[1] + el_j*dx[1] + ((lCoord[1] + 1.0)/2.0)*dx[1];
	}

	/* global shape func derivs */
	gGNix = new double*[2];
	gGNix[0] = new double[el->nNodes];
	gGNix[1] = new double[el->nNodes];
	gGNixx = new double*[3];
	gGNixx[0] = new double[el->nNodes];
	gGNixx[1] = new double[el->nNodes];
	gGNixx[2] = new double[el->nNodes];

	elNodes = new int[el->nNodes];
}

Mesh::~Mesh() {
	for( int vert_i = 0; vert_i < nVertsTotal; vert_i++ ) {
		delete[] verts[vert_i];
	}
	delete[] verts;

	delete[] elNodes;

	delete[] gGNix[0];
	delete[] gGNix[1];
	delete[] gGNix;

	delete[] gGNixx[0];
	delete[] gGNixx[1];
	delete[] gGNixx[2];
	delete[] gGNixx;

	delete el;
}

#define EPSILON 1.0e-12
#define NUMAPPROX( var, val ) ( var >= val - EPSILON && var <= val + EPSILON )

void Mesh::GetEl( double* coord, int* el_i ) {
        int topo[2];
	double  out, frac, integer;

        if( NUMAPPROX( coord[0] - max[0], 0.0 ) ) {
                topo[0] = nEls[0] - 1;
        }
        else {
                out = (coord[0] - min[0])/dx[0];
                frac = modf( out, &integer );
                topo[0] = (int)integer;
                if( topo[0] > 0 && NUMAPPROX( frac, 0.0 ) ) {
                        topo[0]--;
                }
        }
        if( NUMAPPROX( coord[1] - max[1], 0.0 ) ) {
                topo[1] = nEls[1] - 1;
        }
        else {
                out = (coord[1] - min[1])/dx[1];
                frac = modf( out, &integer );
                topo[1] = (int)integer;
                if( topo[1] > 0 && NUMAPPROX( frac, 0.0 ) ) {
                        topo[1]--;
                }
        }
        *el_i = topo[1]*nEls[0] + topo[0];
}

int* Mesh::ElNodes( int el_i ) {
        int topo[2], ii, jj, node_i = 0, n0;
	
	topo[0] = N*(el_i%nEls[0]);
	topo[1] = N*(el_i/nEls[0]);
	n0 = topo[1]*nVerts[0] + topo[0];

	for( jj = 0; jj < N+1; jj++ ) {
		for( ii = 0; ii < N+1; ii++ ) {
			elNodes[node_i++] = n0 + jj*nVerts[0] + ii;
		}
	}

	return elNodes;
}

void Mesh::IndexToTopo( int index, int* topo ) {
	topo[0] = index%nVerts[0];
	topo[1] = index/nVerts[0];
}

void Mesh::TopoToIndex( int* topo, int* index ) {
	*index = topo[1]*nVerts[0] + topo[0];
}

void Mesh::GlobalToLocal( double* gCoord, int* el_i, double* lCoord ) {
	int dim_i;
	GetEl( gCoord, el_i );
	lCoord[0] = (gCoord[0] - min[0] - (*el_i%nEls[0])*dx[0])*(2.0/dx[0]) - 1.0;
	lCoord[1] = (gCoord[1] - min[1] - (*el_i/nEls[0])*dx[1])*(2.0/dx[1]) - 1.0;
	for( dim_i = 0; dim_i < 2; dim_i++ ) {
		if( lCoord[dim_i] < -1.0 ) {
			lCoord[dim_i] = -1;
		}
		if( lCoord[dim_i] > +1.0 ) {
			lCoord[dim_i] = +1;
		}
	}
}

void Mesh::LocalToGlobal( double* lCoord, int el_i, double* gCoord ) {
	gCoord[0] = min[0] + (el_i%nEls[0])*dx[0] + 0.5*(lCoord[0] + 1.0)*dx[0];
	gCoord[1] = min[1] + (el_i/nEls[0])*dx[1] + 0.5*(lCoord[1] + 1.0)*dx[1];
}

double Mesh::DetJac( int el_i, int pt_i ) {
	double J[2][2];
	double** GNix = el->ShapeFuncDerivs( pt_i );

	J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
	ElNodes( el_i );

	for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
		J[0][0] += GNix[0][node_i]*verts[elNodes[node_i]][0];
		J[0][1] += GNix[0][node_i]*verts[elNodes[node_i]][1];
		J[1][0] += GNix[1][node_i]*verts[elNodes[node_i]][0];
		J[1][1] += GNix[1][node_i]*verts[elNodes[node_i]][1];
	}

	return J[0][0]*J[1][1] - J[0][1]*J[1][0];
}

double Mesh::DetJacAtCoord( int el_i, double* coord ) {
	double J[2][2];
	double** GNix = el->ShapeFuncDerivsAtCoord( coord );

	J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
	ElNodes( el_i );

	for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
		J[0][0] += GNix[0][node_i]*verts[elNodes[node_i]][0];
		J[0][1] += GNix[0][node_i]*verts[elNodes[node_i]][1];
		J[1][0] += GNix[1][node_i]*verts[elNodes[node_i]][0];
		J[1][1] += GNix[1][node_i]*verts[elNodes[node_i]][1];
	}

	return J[0][0]*J[1][1] - J[0][1]*J[1][0];
}

double** Mesh::ShapeFuncDerivs( int el_i, int pt_i, double* detJac ) {
	double J[2][2], temp;
	double** GNix = el->ShapeFuncDerivs( pt_i );

	ElNodes( el_i );
	
	J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
	for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
		J[0][0] += GNix[0][node_i]*verts[elNodes[node_i]][0];
		J[0][1] += GNix[0][node_i]*verts[elNodes[node_i]][1];
		J[1][0] += GNix[1][node_i]*verts[elNodes[node_i]][0];
		J[1][1] += GNix[1][node_i]*verts[elNodes[node_i]][1];
	}
	*detJac = J[0][0]*J[1][1] - J[0][1]*J[1][0];

	/* invert J */
	temp = J[0][0];
	J[0][0] = J[1][1]/(*detJac);
	J[1][1] = temp/(*detJac);
	J[0][1] = -J[0][1]/(*detJac);
	J[1][0] = -J[1][0]/(*detJac);
	
	for( int dim_i = 0; dim_i < 2; dim_i++ ) {
		for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
			temp = 0.0;
			for( int dim_j = 0; dim_j < 2; dim_j++ ) {
				temp += GNix[dim_j][node_i]*J[dim_i][dim_j];
			}
			gGNix[dim_i][node_i] = temp;
		}
	}

	return gGNix;
}

double** Mesh::ShapeFuncDerivsAtCoord( int el_i, double* coord, double* detJac ) {
	double J[2][2], temp;
	double** GNix = el->ShapeFuncDerivsAtCoord( coord );

	ElNodes( el_i );
	
	J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
	for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
		J[0][0] += GNix[0][node_i]*verts[elNodes[node_i]][0];
		J[0][1] += GNix[0][node_i]*verts[elNodes[node_i]][1];
		J[1][0] += GNix[1][node_i]*verts[elNodes[node_i]][0];
		J[1][1] += GNix[1][node_i]*verts[elNodes[node_i]][1];
	}
	*detJac = J[0][0]*J[1][1] - J[0][1]*J[1][0];

	/* invert J */
	temp = J[0][0];
	J[0][0] = J[1][1]/(*detJac);
	J[1][1] = temp/(*detJac);
	J[0][1] = -J[0][1]/(*detJac);
	J[1][0] = -J[1][0]/(*detJac);
	
	for( int dim_i = 0; dim_i < 2; dim_i++ ) {
		for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
			temp = 0.0;
			for( int dim_j = 0; dim_j < 2; dim_j++ ) {
				temp += GNix[dim_j][node_i]*J[dim_i][dim_j];
			}
			gGNix[dim_i][node_i] = temp;
		}
	}

	return gGNix;
}

double** Mesh::ShapeFuncSecondDerivs( int el_i, int pt_i, double* detJac ) {
/*
	double 		**GNix, **GNixx, j2Inv, j3Inv, ja, jb;
	double 		fr, frr, fs, fss, gr, grr, gs, gss, frs, grs;
	Legendre*	l 	= (Legendre*)el;

	ElNodes( el_i );
	ShapeFuncDerivs( el_i, pt_i, detJac );
	GNix  = l->ShapeFuncDerivs( pt_i );
	GNixx = l->ShapeFunc2ndDerivs( pt_i );

	fr = frr = fs = fss = gr = grr = gs = gss = frs = grs = 0.0;
	for( int node_i = 0; node_i < l->nNodes; node_i++ ) {
		fr  += GNix[0][node_i]*verts[elNodes[node_i]][0];
		fs  += GNix[1][node_i]*verts[elNodes[node_i]][0];
		gr  += GNix[0][node_i]*verts[elNodes[node_i]][1];
		gs  += GNix[1][node_i]*verts[elNodes[node_i]][1];
		frr += GNixx[0][node_i]*verts[elNodes[node_i]][0];
		fss += GNixx[1][node_i]*verts[elNodes[node_i]][0];
		frs += GNixx[2][node_i]*verts[elNodes[node_i]][0];
		grr += GNixx[0][node_i]*verts[elNodes[node_i]][1];
		gss += GNixx[1][node_i]*verts[elNodes[node_i]][1];
		grs += GNixx[2][node_i]*verts[elNodes[node_i]][1];
	}

	j2Inv = 1.0/((*detJac)*(*detJac));
	j3Inv = 1.0/((*detJac)*(*detJac)*(*detJac));
	ja = frr*gs + fr*grs - frs*gr - fs*grr;
	jb = frs*gs + fr*gss - fss*gr - fs*grs;
	for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
		gGNixx[0][node_i] = gs*gs*GNixx[0][node_i] - 2.0*gr*gs*GNixx[2][node_i] + gr*gr*GNixx[1][node_i] -
				    gGNix[1][node_i]*( gs*gs*grr - 2.0*gr*gs*grs + gr*gr*gss ) - 
				    gGNix[0][node_i]*( gs*gs*frr - 2.0*gr*gs*fss + gr*gr*fss );
		gGNixx[0][node_i] *= j2Inv;

		gGNixx[1][node_i] = fs*fs*GNixx[0][node_i] - 2.0*fr*fs*GNixx[2][node_i] + fr*fr*GNixx[1][node_i] -
				    gGNix[1][node_i]*( fs*fs*grr - 2.0*fr*fs*grs + fr*fr*gss ) - 
				    gGNix[0][node_i]*( fs*fs*frr - 2.0*fr*fs*fss + fr*fr*fss );
		gGNixx[1][node_i] *= j2Inv;

		gGNixx[2][node_i] = j2Inv*( (fr*gs + fs*gr)*GNixx[2][node_i] - fr*gr*GNixx[1][node_i] - fs*gs*GNixx[0][node_i] ) + 
				    ( j2Inv*(fr*gss - fs*grs) + j3Inv*(fs*gs*ja - fr*gs*jb) )*GNix[0][node_i] +
				    ( j2Inv*(fs*grr - fr*grs) + j3Inv*(fr*gr*jb - fs*gr*ja) )*GNix[1][node_i];
	}
	return gGNixx;
*/
	double **GNixx, a, b;
	Legendre* l = (Legendre*)el;

	a = 2.0/dx[0];
	b = 2.0/dx[1];
	*detJac = 1.0/a/b;
	GNixx = l->ShapeFuncSecondDerivs( pt_i );
	for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
		gGNixx[0][node_i] = a*a*GNixx[0][node_i];
		gGNixx[1][node_i] = b*b*GNixx[1][node_i];
		gGNixx[2][node_i] = 0.0;
	}
	return gGNixx;
}

double** Mesh::ShapeFuncSecondDerivsAtCoord( int el_i, double* coord, double* detJac ) {
	double** 	GNixx, a, b;
	Legendre*	l	= (Legendre*)el;

	a = 2.0/dx[0];
	b = 2.0/dx[1];
	*detJac = 1.0/a/b;
	GNixx = l->ShapeFuncSecondDerivsAtCoord( coord );
	for( int node_i = 0; node_i < el->nNodes; node_i++ ) {
		gGNixx[0][node_i] = a*a*GNixx[0][node_i];
		gGNixx[1][node_i] = b*b*GNixx[1][node_i];
		gGNixx[2][node_i] = 0.0;
	}
	return gGNixx;
}

//#define HDF5_OLD

void Mesh::Save() {
	hid_t 		file, attrib_id, group_id, attribData_id, fileSpace, fileData, fileSpace2, fileData2, memSpace;
	hsize_t 	a_dims, start[2], count[2], size[2];
	int		attribData;
	char 		filename[40];
	int 		nLinearEls[2];
	int 		nDims			= 2;
	int		nodesPerLinearEl	= 4;
	int		node_i, el_i;

	sprintf( filename, "%s.h5", name.c_str() );
	file = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

	/* mesh dimensionality */
	a_dims = 1;
	attribData = 2; /* mesh dimensionality */
	attribData_id = H5Screate_simple( 1, &a_dims, NULL );
#ifdef HDF5_OLD
	group_id = H5Gopen( file, "/" );
	attrib_id = H5Acreate( group_id, "dimensions", H5T_STD_I32BE, attribData_id, H5P_DEFAULT );
#else
	group_id = H5Gopen( file, "/", H5P_DEFAULT );
	attrib_id = H5Acreate( group_id, "dimensions", H5T_STD_I32BE, attribData_id, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Awrite( attrib_id, H5T_NATIVE_INT, &attribData );
	H5Aclose( attrib_id );
	H5Gclose( group_id );
	H5Sclose( attribData_id );
	
	/* mesh resolution */
	nLinearEls[0] = nVerts[0] - 1;
	nLinearEls[1] = nVerts[1] - 1;
	a_dims = 2;
	attribData_id = H5Screate_simple( 1, &a_dims, NULL );
#ifdef HDF5_OLD
	group_id = H5Gopen( file, "/" );
	attrib_id = H5Acreate( group_id, "mesh resolution", H5T_STD_I32BE, attribData_id, H5P_DEFAULT );
#else
	group_id = H5Gopen( file, "/", H5P_DEFAULT );
	attrib_id = H5Acreate( group_id, "mesh resolution", H5T_STD_I32BE, attribData_id, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Awrite( attrib_id, H5T_NATIVE_INT, nLinearEls );
	H5Aclose( attrib_id );
	H5Gclose( group_id );
	H5Sclose( attribData_id );

	/* max and min coords of mesh */
	count[0] = (hsize_t)nDims;
	fileSpace = H5Screate_simple( 1, count, NULL );
#ifdef HDF5_OLD
	fileData = H5Dcreate( file, "/min", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
#else
	fileData = H5Dcreate( file, "/min", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Dwrite( fileData, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, min );
	H5Dclose( fileData );
	H5Sclose( fileSpace );
	fileSpace = H5Screate_simple( 1, count, NULL );
#ifdef HDF5_OLD
	fileData = H5Dcreate( file, "/max", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
#else
	fileData = H5Dcreate( file, "/max", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Dwrite( fileData, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, max );
	H5Dclose( fileData );
	H5Sclose( fileSpace );

	/* write the vertices */
 	size[0] = (hsize_t)nVertsTotal;
	size[1] = (hsize_t)nDims;
	fileSpace = H5Screate_simple( 2, size, NULL );
#ifdef HDF5_OLD
	fileData = H5Dcreate( file, "/vertices", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
#else
	fileData = H5Dcreate( file, "/vertices", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif
	size[0] = (hsize_t)nLinearEls[0]*nLinearEls[1]; /* use linear quads for writing file */
	size[1] = (hsize_t)nodesPerLinearEl;
	fileSpace2 = H5Screate_simple( 2, size, NULL );
#ifdef HDF5_OLD
	fileData2 = H5Dcreate( file, "/connectivity", H5T_NATIVE_INT, fileSpace2, H5P_DEFAULT );
#else
	fileData2 = H5Dcreate( file, "/connectivity", H5T_NATIVE_INT, fileSpace2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif

	count[0] = 1;
	count[1] = nDims;
	memSpace = H5Screate_simple( 2, count, NULL );
	H5Sselect_all( memSpace );

	for( node_i = 0; node_i < nVertsTotal; node_i++ ) {
		start[1] = 0;
		start[0] = node_i;
		H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Dwrite( fileData, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, verts[node_i] );
	}
	H5Sclose( memSpace );
	H5Dclose( fileData );
	H5Sclose( fileSpace );

	H5Sget_simple_extent_dims( fileSpace2, size, NULL );
	count[0] = 1;
	count[1] = size[1];
	memSpace = H5Screate_simple( 2, count, NULL );
	H5Sselect_all( memSpace );
	
	for( el_i = 0; el_i < nLinearEls[0]*nLinearEls[1]; el_i++ ) {
		int linearElNodes[4], j;
		j = el_i/nLinearEls[0];
		linearElNodes[0] = el_i + j;
		linearElNodes[1] = linearElNodes[0] + 1;
		linearElNodes[2] = linearElNodes[1] + nVerts[0];
		linearElNodes[3] = linearElNodes[2] - 1;	
		start[1] = 0;
		start[0] = el_i;
		H5Sselect_hyperslab( fileSpace2, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Dwrite( fileData2, H5T_NATIVE_INT, memSpace, fileSpace2, H5P_DEFAULT, linearElNodes );
	}

	H5Sclose( memSpace );
	H5Dclose( fileData2 );
	H5Sclose( fileSpace2 );
	H5Fclose( file );
}
