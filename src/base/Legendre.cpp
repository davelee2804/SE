#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Legendre.h"

using namespace std;

#define STEP 1.0e-3

Legendre::Legendre( int _N ) : Element( _N ) {
	int endLoop, i, j;
	double x0, *coord;
	polys = new Jacobi*[N+2];
	polys[0] = new Jacobi( 0, 0.0 );
	polys[0]->c[0] = 1.0;
	polys[1] = new Jacobi( 1, 0.0 );
	polys[1]->c[1] = 1.0;
	for( i = 2; i <= N+1; i++ ) {
		polys[i] = new Jacobi( i, 0.0 );
		polys[i]->Gen( polys[i-1], polys[i-2] );
	}

	abcissa[0] = -1.0;
        /* find the roots of the derivative of the Nth Legendre polynomial - these are the Gauss-Lobatto
         * quadrature points... for an initial guess close to the previous root, the previous root will 
         * be returned, so need to keep iterating through till the root found is the next one along, and
         * check against the newton-raphson solver overshooting. only the first half need be calculated.
         * the rest can be obtained through symmetry */
        endLoop = ( (N+1)%2 == 1 ) ? N/2 : (N + 1)/2;
        for( i = 1; i < endLoop; i++ ) {
                x0 = abcissa[i-1];
                do {
                        x0 += STEP;
                        abcissa[i] = RootFinder( x0, polys[N] );
                } while( AlreadyFound( i, abcissa, abcissa[i] ) || abcissa[i] > 0.0 - STEP );
                abcissa[N-i] = -abcissa[i];
        }
        abcissa[N] = 1.0;
        if( (N+1)%2 == 1 ) {
                abcissa[(N+1)/2] = 0.0;
        }

	SortEntries( abcissa, N+1 );
	weights[0] = 2.0/((N+1)*N);
	for( i = 1; i < N; i++ ) {
		weights[i] = 2.0/((N+1)*N*polys[N]->Eval(abcissa[i])*polys[N]->Eval(abcissa[i]));
	}
	weights[N] = 2.0/((N+1)*N);

	for( j = 0; j <= N; j++ ) {
		for( i = 0; i <= N; i++ ) {
			coord = new double[2];
			coord[0] = abcissa[i];
			coord[1] = abcissa[j];
			quadPts[j*(N+1) + i] = new QuadPoint( coord, weights[i]*weights[j] );

			xi[j*(N+1) + i] = new double[2];
			xi[j*(N+1) + i][0] = abcissa[i];
			xi[j*(N+1) + i][1] = abcissa[j];
		}
	}

	P = new Jacobi( N, 0.0 );
	P->Gen( polys[N-1], polys[N-2] );

	dCij  = new double[(N+1)*(N+1)];
	d2Cij = new double[(N+1)*(N+1)];
	GNixx = new double*[3];
	GNixx[0] = new double[nNodes]; /* Ni,xx */
	GNixx[1] = new double[nNodes]; /* Ni,yy */
	GNixx[2] = new double[nNodes]; /* Ni,xy */

	/* cardinal function derivative matrix,
		i: quadrature point index
		j: basis function index         
	   references:
		Boyd (2000) Chebychev and Fourier Spectral Methods, Second Ed. Dover 
		Gottlieb, Hussaini and Orszag (1984) Theory and Applications of Spectral Methods, 
			Spectral Methods for Partial Differential Equaitons, pp 1-54, Voigt, 
			Gottlieb and Hussaini, eds. SIAM */
	for( i = 0; i <= N; i++ ) {
		for( j = 0; j <= N; j++ ) {
                        if     ( i == j && i == 0 ) { dCij[i*(N+1)+j] = -0.25*N*(N + 1); }
                        else if( i == j && i == N ) { dCij[i*(N+1)+j] = +0.25*N*(N + 1); }
                        else if( i == j )           { dCij[i*(N+1)+j] = 0.0; }
                        else { 
                                dCij[i*(N+1)+j] = P->Eval( abcissa[i] )/( P->Eval( abcissa[j] )*( abcissa[i] - abcissa[j] ) ); 
                        }
		}
	}
	/* cardinal function second derivative matrix,
		i: quadrature point index
		j: basis function index         
	   references:
		Gottlieb, Hussaini and Orszag (1984) Theory and Applications of Spectral Methods, 
			Spectral Methods for Partial Differential Equaitons, pp 1-54, Voigt, 
			Gottlieb and Hussaini, eds. SIAM 
		Solomonoff and Turkel (1989) Global Properties of Pseudospectral Methods,
			J. Comp. Phys. 81, 239-276 */
	for( i = 0; i <= N; i++ ) {
		for( j = 0; j <= N; j++ ) {
			if( i == j ) {
				d2Cij[i*(N+1)+j] = dCij[i*(N+1)+j]*dCij[i*(N+1)+j];
				for( int k = 0; k <= N; k++ ) {
					if( k != i ) {
						d2Cij[i*(N+1)+j] -= 1.0/(abcissa[k] - abcissa[i])/(abcissa[k] - abcissa[i]);
					}
				}
			}
			else {
				d2Cij[i*(N+1)+j] = 0.0;
				for( int k = 0; k <= N; k++ ) {
					if( k != i && k != j ) {
						d2Cij[i*(N+1)+j] += 2.0*dCij[i*(N+1)+j]/(abcissa[i] - abcissa[k]);
					}
				}
			}
		}
	}
}

Legendre::~Legendre() {
	for( int i = 0; i <= N+1; i++ ) { delete polys[i]; }
	delete[] polys;
	delete[] dCij;
	delete[] d2Cij;
	delete[] GNixx[0];
	delete[] GNixx[1];
	delete[] GNixx[2];
	delete[] GNixx;
}

double* Legendre::ShapeFuncs( int pt ) {
	for( int i = 0; i < nNodes; i++ ) { Ni[i] = 0.0; }
	Ni[pt] = 1.0;

	return Ni;
}

double* Legendre::ShapeFuncsAtCoord( double* coord ) {
	int x_i, y_i, x_j, y_j;
	double Nx, Ny;

	for( int pt_j = 0; pt_j < nNodes; pt_j++ ) {
		x_j = pt_j%(N+1);
		y_j = pt_j/(N+1);

		if( IsRoot( coord[0], &x_i ) ) {
			Nx = ( x_i == x_j ) ? 1.0 : 0.0;
		}
		else {
			Nx = -((1.0 - coord[0]*coord[0])/(N*(N + 1)*P->Eval(abcissa[x_j])*(coord[0] - abcissa[x_j])))*P->EvalDeriv(coord[0]);
		}

		if( IsRoot( coord[1], &y_i ) ) {
			Ny = ( y_i == y_j ) ? 1.0 : 0.0;
		}
		else {
			Ny = -((1.0 - coord[1]*coord[1])/(N*(N + 1)*P->Eval(abcissa[y_j])*(coord[1] - abcissa[y_j])))*P->EvalDeriv(coord[1]);
		}
		Ni[pt_j] = Nx*Ny;
	}
	return Ni;
}

double** Legendre::ShapeFuncDerivs( int pt ) {
        int x_i, y_i, x_j, y_j;

	x_i = pt%(N+1);
	y_i = pt/(N+1);

        for( int pt_j = 0; pt_j < nNodes; pt_j++ ) {
                x_j = pt_j%(N+1);
                y_j = pt_j/(N+1);
		GNix[0][pt_j] = ( y_i == y_j ) ? dCij[x_i*(N+1)+x_j] : 0.0;
		GNix[1][pt_j] = ( x_i == x_j ) ? dCij[y_i*(N+1)+y_j] : 0.0;
        }
	return GNix;
}

// TODO: check this - even though is not currently in use...
double** Legendre::ShapeFuncDerivsAtCoord( double* coord ) {
	double Nx, Ny, GNx, GNy, a, b, c;
	int x_i, y_i, x_j, y_j;

	for( int pt_j = 0; pt_j < nNodes; pt_j++ ) {
		x_j = pt_j%(N+1);
		y_j = pt_j/(N+1);

		if( IsRoot( coord[0], &x_i ) ) {
			Nx = ( x_i == x_j ) ? 1.0 : 0.0;
			GNx = dCij[x_i*(N+1)+x_j];
		}
		else {
			Nx = -((1.0 - coord[0]*coord[0])/(N*(N + 1)*P->Eval(abcissa[x_j])*(coord[0] - abcissa[x_j])))*P->EvalDeriv(coord[0]);
	        	a = coord[0]*coord[0] + 1.0 - 2.0*coord[0]*abcissa[x_j];
                	b = (1.0 - coord[0]*coord[0])*(coord[0] - abcissa[x_j]);
                	c = N*(N + 1)*(coord[0] - abcissa[x_j])*(coord[0] - abcissa[x_j])*P->Eval(abcissa[x_j]);
                	GNx = (a*P->EvalDeriv(coord[0]) - b*P->Eval2ndDeriv(coord[0]))/c;
		}

		if( IsRoot( coord[1], &y_i ) ) {
			Ny = ( y_i == y_j ) ? 1.0 : 0.0;
			GNy = dCij[y_i*(N+1)+y_j];
		}
		else {
			Ny = -((1.0 - coord[1]*coord[1])/(N*(N + 1)*P->Eval(abcissa[y_j])*(coord[1] - abcissa[y_j])))*P->EvalDeriv(coord[1]);
	        	a = coord[1]*coord[1] + 1.0 - 2.0*coord[1]*abcissa[y_j];
                	b = (1.0 - coord[1]*coord[1])*(coord[1] - abcissa[y_j]);
                	c = N*(N + 1)*(coord[1] - abcissa[y_j])*(coord[1] - abcissa[y_j])*P->Eval(abcissa[y_j]);
                	GNy = (a*P->EvalDeriv(coord[1]) - b*P->Eval2ndDeriv(coord[1]))/c;
		}

		GNix[0][pt_j] = GNx*Ny;
		GNix[1][pt_j] = Nx*GNy;
	}
	return GNix;
}

double** Legendre::ShapeFuncSecondDerivsAtCoord( double* x ) {
	double 	Nx, Ny, GNx, GNy, GNxx, GNyy;
	double	xInv, phi, phij, dphi, alpha;
	int 	x_i, y_i, x_j, y_j;

	for( int pt_j = 0; pt_j < nNodes; pt_j++ ) {
		x_j = pt_j%(N+1);
		y_j = pt_j/(N+1);

		if( IsRoot( x[0], &x_i ) ) {
			Nx   = ( x_i == x_j ) ? 1.0 : 0.0;
			GNx  = dCij[x_i*(N+1)+x_j];
			GNxx = d2Cij[x_i*(N+1)+x_j];
		}
		else {
			xInv  = 1.0/(x[0] - abcissa[x_j]);
			phi   = P->Eval( x[0] );
			phij  = P->Eval( abcissa[x_j] );
			dphi  = P->EvalDeriv( x[0] );
			alpha = 1.0/(N*(N+1)*phij);

			Nx   = (x[0]*x[0] - 1.0)*dphi*xInv*alpha;
			GNx  = phi*xInv/phij + (1.0 - x[0]*x[0])*dphi*xInv*xInv*alpha;
			GNxx = dphi*xInv/phij - 2.0*phi*xInv*xInv/phij - 2.0*(1.0 - x[0]*x[0])*dphi*xInv*xInv*xInv*alpha;
		}
		if( IsRoot( x[1], &y_i ) ) {
			Ny   = ( y_i == y_j ) ? 1.0 : 0.0;
			GNy  = dCij[y_i*(N+1)+y_j];
			GNyy = d2Cij[y_i*(N+1)+y_j];
		}
		else {
			xInv  = 1.0/(x[1] - abcissa[y_j]);
			phi   = P->Eval( x[1] );
			phij  = P->Eval( abcissa[y_j] );
			dphi  = P->EvalDeriv( x[1] );
			alpha = 1.0/(N*(N+1)*phij);

			Ny   = (x[1]*x[1] - 1.0)*dphi*xInv*alpha;
			GNy  = phi*xInv/phij + (1.0 - x[1]*x[1])*dphi*xInv*xInv*alpha;
			GNyy = dphi*xInv/phij - 2.0*phi*xInv*xInv/phij - 2.0*(1.0 - x[1]*x[1])*dphi*xInv*xInv*xInv*alpha;
		}

		GNixx[0][pt_j] = GNxx*Ny;
		GNixx[1][pt_j] = Nx*GNyy;
		GNixx[2][pt_j] = GNx*GNy;
	}

	return GNixx;
}

double** Legendre::ShapeFuncSecondDerivs( int pt ) {
	int x_i, y_i, x_j, y_j;

	x_i = pt%(N+1);
	y_i = pt/(N+1);

	for( int pt_j = 0; pt_j < nNodes; pt_j++ ) {
		x_j = pt_j%(N+1);
		y_j = pt_j/(N+1);
		GNixx[0][pt_j] = ( y_i == y_j ) ? d2Cij[x_i*(N+1)+x_j] : 0.0;
		GNixx[1][pt_j] = ( x_i == x_j ) ? d2Cij[y_i*(N+1)+y_j] : 0.0;
		GNixx[2][pt_j] = dCij[x_i*(N+1)+x_j]*dCij[y_i*(N+1)+y_j];
	}
	return GNixx;
}

#define EPS 1.0e-8

bool Legendre::IsRoot( double x, int* i ) {
        int k;

        for( k = 0; k <= N; k++ ) {
                if( fabs( x - abcissa[k] ) < EPS ) {
                        *i = k;
                        return true;
                }
        }

        *i = -1;
        return false;
}

Jacobi** Legendre::ModalBasis() {
        Jacobi**        Pn;
        Jacobi**        modalBasis;

        Pn = new Jacobi*[N+1];
        Pn[0] = new Jacobi( 0, 1 );
        Pn[0]->c[0] = 1.0;
        Pn[1] = new Jacobi( 1, 1 );
        Pn[1]->c[1] = 2.0;
        for( int i = 2; i <= N; i++ ) {
                Pn[i] = new Jacobi( i, 1 );
                Pn[i]->Gen( Pn[i-1], Pn[i-2] );
        }

        modalBasis = new Jacobi*[N+1];
        modalBasis[0] = new Jacobi( 1, 1 );
        modalBasis[0]->c[0] = +0.5;
        modalBasis[0]->c[1] = +0.5;
        modalBasis[1] = new Jacobi( 1, 1 );
        modalBasis[1]->c[0] = +0.5;
        modalBasis[1]->c[1] = -0.5;
        for( int i = 2; i <= N; i++ ) {
                modalBasis[i] = new Jacobi( i, 1 );
                for( int j = 0; j < i-1; j++ ) {
                        modalBasis[i]->c[j] += 0.25*Pn[i-2]->c[j];
                        modalBasis[i]->c[j+2] -= 0.25*Pn[i-2]->c[j];
                }
        }
        for( int i = 0; i <= N; i++ ) { delete Pn[i]; }
        delete[] Pn;

        return modalBasis;
}

double* Legendre::ModalToNodalTransformMatrix() {
	double*		Bij	= new double[nNodes*nNodes];
	Jacobi**	mBasis	= ModalBasis();
	int		xi, yi, xj, yj;

        for( int i = 0; i < nNodes; i++ ) {
                xi = i%(N+1);
                yi = i/(N+1);
                for( int j = 0; j < nNodes; j++ ) {
                        xj = j%(N+1);
                        yj = j/(N+1);
                        Bij[i*nNodes+j] = mBasis[xj]->EvalTotal( abcissa[xi] )*mBasis[yj]->EvalTotal( abcissa[yi] );
                }
        }

	for( int i = 0; i <= N; i++ ) { delete mBasis[i]; }
	delete[] mBasis;

	return Bij;
}
