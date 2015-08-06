#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Legendre.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Operator.h"
#include "SVV.h"

using namespace std;
using std::string;

#define I(p,q,N) \
	p*N + q

void GetRowCol( int i, int* p, int* q, int N ) { *p = i/N; *q = i%N; }

/* spectral vanishing viscosity (SVV) for filtering out high frequency modes.
 * configure the petsc solver as cg/jacobi with rtol = 1.0e-20 or higher.
 * references: 
 *	Kirby, R. M. & S. J. Sherwin (2006) "Stabilisation of spectral/hp 
 *	element methods through spectral vanishing viscosity: Application
 *	to fluid mechanics modelling" Comput. Methods Appl. Mech. Engrg. 
 *	195, 3128 - 3144
 *      Xu, C. & Pasquetti R. (2004) "Stabilized spectral element computations
 *      of high Reynolds number incompressible flows" J. Comput. Phys. 196,
 *      680 - 704
 *      Blackburn, H. M. & S. Schmidt (2003) "Spectral element filtering
 *      techniques for large eddy simulation with dynamic estimation" 
 *	J. Comput. Phys. 186, 610 - 629
 *	Karamanos, G. S. and G. E. Karniadakis (2000) "A Spectral Vanishing
 *	Viscosity Method for Large-Eddy Simulations" J. Comput. Phys. 163,
 *	22 - 50
 */

SVV::SVV( string _name, Field* _rowField, Field* _colField, double _constant, int _mN ) : Operator( _name, _rowField, _colField, _constant ) {
	Legendre* 	el 	= (Legendre*)_rowField->mesh->el;
	int 		nNodes 	= el->nNodes;
	int 		nPoints = el->nPoints;

	mN = _mN;

	Q			= new double[nNodes*nNodes];
	W			= new double[nNodes*nNodes];
	Lxx			= new double[nNodes*nNodes];
	Lyy			= new double[nNodes*nNodes];
	B			= new double[nPoints*nNodes];
	Bt			= new double[nNodes*nPoints];
	BtW			= new double[nNodes*nNodes];
	BtWB			= new double[nPoints*nPoints];
	BtWBinv			= new double[nPoints*nPoints];
	M			= new double[nNodes*nNodes];
	Minv			= new double[nNodes*nNodes];
	Dx			= new double[nPoints*nNodes];
	Dy			= new double[nPoints*nNodes];
	Dxt			= new double[nNodes*nPoints];
	Dyt			= new double[nNodes*nPoints];
	MinvQ			= new double[nNodes*nNodes];
	MinvQM			= new double[nNodes*nNodes];
	MinvQMDx		= new double[nNodes*nNodes];
	MinvQMDy		= new double[nNodes*nNodes];
	MinvQMDxt		= new double[nNodes*nNodes];
	MinvQMDyt		= new double[nNodes*nNodes];
	MinvQMDxtW		= new double[nNodes*nNodes];
	MinvQMDytW		= new double[nNodes*nNodes];

	mBasis = el->ModalBasis();

	if( rowField != colField ) { 
		cerr << "ERROR: attempting to apply SVV operator where row field != column field.\n";
		exit(1);
	}
}

SVV::~SVV() {
	for( int i = 0; i <= rowField->mesh->el->N; i++ ) {
		delete mBasis[i];
	}
	delete[] mBasis;

	delete[] Q;
	delete[] W;
	delete[] Lxx;
	delete[] Lyy;
	delete[] B;
	delete[] Bt;
	delete[] BtW;
	delete[] BtWB;
	delete[] BtWBinv;
	delete[] M;
	delete[] Minv;
	delete[] Dx;
	delete[] Dy;
	delete[] Dxt;
	delete[] Dyt;
	delete[] MinvQ;
	delete[] MinvQM;
	delete[] MinvQMDx;
	delete[] MinvQMDy;
	delete[] MinvQMDxt;
	delete[] MinvQMDyt;
	delete[] MinvQMDxtW;
	delete[] MinvQMDytW;
}

void SVV::AssembleElement( int el_i, double* A ) {
	Mesh*		mesh		= rowField->mesh;
	Element*	el		= mesh->el;
	double		**GNix, detJac, weight, qExpX, qExpY;
	int		ix, iy, jx, jy, row, col;
	int		nDofs = rowField->nDofs, nNodes = el->nNodes, nEntries = nDofs*nNodes;

	for( int i = 0; i < nNodes; i++ ) {
		GetRowCol( i, &ix, &iy, el->N+1 );
		weight = el->quadPts[i]->weight;
		GNix   = el->ShapeFuncDerivs( i );
		detJac = mesh->DetJac( el_i, i );
		for( int j = 0; j < nNodes; j++ ) {
			GetRowCol( j, &jx, &jy, el->N+1 );
			B[I(i,j,nNodes)] = mBasis[jx]->EvalTotal( el->abcissa[ix] )*mBasis[jy]->EvalTotal( el->abcissa[iy] );
			W[I(i,j,nNodes)] = Q[I(i,j,nNodes)] = 0.0;
			if( i == j ) {
				W[I(i,j,nNodes)] = detJac*weight;
				if( jx > mN && jy > mN ) {
					qExpX = -1.0*( (double)(el->N - jx)*(el->N - jx) )/( (mN - jx)*(mN - jx) );
					qExpY = -1.0*( (double)(el->N - jy)*(el->N - jy) )/( (mN - jy)*(mN - jy) );
					Q[I(i,j,nNodes)] = sqrt( exp( qExpX + qExpY ) );
				}
			}
			Dx[I(i,j,nNodes)] = GNix[0][j];
			Dy[I(i,j,nNodes)] = GNix[1][j];
		}
	}
	Tran( Dx, Dxt, nNodes );
	Tran( Dy, Dyt, nNodes );
	Tran( B, Bt, nNodes );
	Mult( Bt, W, BtW, nNodes );
	Mult( BtW, B, BtWB, nNodes );
	Inv( BtWB, BtWBinv, nNodes );
	Mult( BtWBinv, BtW, M, nNodes );
	Inv( M, Minv, nNodes );
	Mult( Minv, Q, MinvQ, nNodes );
	Mult( MinvQ, M, MinvQM, nNodes );
	Mult( MinvQM, Dx, MinvQMDx, nNodes );
	Mult( MinvQM, Dy, MinvQMDy, nNodes );
	Tran( MinvQMDx, MinvQMDxt, nNodes );
	Tran( MinvQMDy, MinvQMDyt, nNodes );
	Mult( MinvQMDxt, W, MinvQMDxtW, nNodes );
	Mult( MinvQMDyt, W, MinvQMDytW, nNodes );
	Mult( MinvQMDxtW, MinvQMDx, Lxx, nNodes );
	Mult( MinvQMDxtW, MinvQMDx, Lyy, nNodes );

	for( int row_i = 0; row_i < nNodes; row_i++ ) {
		for( int col_j = 0; col_j < nNodes; col_j++ ) {
			row = row_i*nDofs;
			col = col_j*nDofs;
			for( int dof = 0; dof < nDofs; dof++ ) {
				A[(row+dof)*nEntries + (col+dof)] += constant*( Lxx[I(row_i,col_j,nNodes)] + Lyy[I(row_i,col_j,nNodes)] );
			}
		}
	}
}
