#include <cmath>
#include <cstdlib>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"

using namespace std;

Element::Element( int _N ) {
	N = _N;

	abcissa = new double[N+1];
	weights = new double[N+1];

	nNodes = (N+1)*(N+1);
	nPoints = (N+1)*(N+1);

	quadPts = new QuadPoint*[nPoints];
	xi = new double*[nNodes];

	Ni = new double[nNodes];
	GNix = new double*[2];
	GNix[0] = new double[nNodes];
	GNix[1] = new double[nNodes];
}

Element::~Element() {
	int pt_i;

	delete[] weights;
	delete[] abcissa;

	delete[] Ni;
	delete[] GNix[0];
	delete[] GNix[1];
	delete[] GNix;

	for( pt_i = 0; pt_i < nPoints; pt_i++ ) {
		delete quadPts[pt_i];
		delete[] xi[pt_i];
	}
	delete[] quadPts;
	delete[] xi;

	delete P;
}

/********************************************************************************************************************/

#define TOLERANCE 1.0E-12

/* newton-raphson method to find the root of a scalar function, given the function and its derivative */
double RootFinder( double x0, Jacobi* P ) {
        double x = x0, xPrevious;

        do {
                xPrevious = x;
                x = xPrevious - P->EvalDeriv(xPrevious)/P->Eval2ndDeriv(xPrevious);
                if( fabs(x - xPrevious) < TOLERANCE ) {
                        return x;
                }
        } while( true );
}

bool AlreadyFound( int n, double* roots, double root ) {
        int i;

        for( i = 0; i < n; i++ ) {
                if( fabs(roots[i] - root) < 1.0e-6 ) {
                        return true;
                }
        }
        return false;
}

void SortEntries( double* vals, int n ) {
	int i, j;
	bool* done = new bool[n];
	double* temp = new double[n];
	double minVal;
	int minInd = 999999999;

	for( i = 0; i < n; i++ ) {
		done[i] = false;
	}

	for( i = 0; i < n; i++ ) {
		minVal = 1.0e+99;
		for( j = 0; j < n; j++ ) {
			if( !done[j] && vals[j] < minVal ) {
				minVal = vals[j];
				minInd = j;
			}
		}
		temp[i] = minVal;
		done[minInd] = true;
	}

	for( i = 0; i < n; i++ ) {
		vals[i] = temp[i];
	}

	delete[] temp;
	delete[] done;
}
/********************************************************************************************************************/

#define STEP 1.0e-3
#define TOL 1.0e-12

Lagrange::Lagrange( int _N ) : Element( _N ) {
	double 		dx 		= 2.0/N;
	bool*		done		= new bool[N+1];
	double*		temp		= new double[N+1];
	int		minInd, endLoop;
	double		minVal, x0, xCurr, xPrev;
	Jacobi**	polys 		= new Jacobi*[N+2];
	double*		coord;

        polys[0] = new Jacobi( 0, 0.0 );
        polys[0]->c[0] = 1.0;
        polys[1] = new Jacobi( 1, 0.0 );
        polys[1]->c[1] = 1.0;
        for( int i = 2; i <= N+1; i++ ) {
                polys[i] = new Jacobi( i, 0.0 );
                polys[i]->Gen( polys[i-1], polys[i-2] );
        }

	P = new Jacobi( N, 0.0 );
	P->Gen( polys[N-1], polys[N-2] );

        /* calculate the gaussian quadrature weights and abcissa */
        endLoop = (N+1)%2 == 1 ? N/2 : (N+1)/2;
        for( int pt_i = 0; pt_i < endLoop; pt_i++ ) {
                x0 = (pt_i == 0) ? -1.0 : abcissa[pt_i-1];
                do {
                        x0 += STEP;
                        xCurr = x0;
                        do {
                                xPrev = xCurr;
                                xCurr = xPrev - P->Eval( xPrev )/P->EvalDeriv( xPrev );
                        } while( fabs( xCurr - xPrev ) > TOL );
                } while( AlreadyFound( pt_i, abcissa, xCurr ) || xCurr > 0.0 - STEP );
                abcissa[pt_i] = xCurr;
                abcissa[N-pt_i] = -abcissa[pt_i];
        }
        if( (N+1)%2 == 1 ) { abcissa[(N+1)/2] = 0.0; }
        minInd = 999999999;
        for( int pt_i = 0; pt_i < N+1; pt_i++ ) { done[pt_i] = false; }
        for( int pt_i = 0; pt_i < N+1; pt_i++ ) {
                minVal = 1.0e+99;
                for( int pt_j = 0; pt_j < N+1; pt_j++ ) {
                        if( !done[pt_j] && abcissa[pt_j] < minVal ) {
                                minVal = abcissa[pt_j];
                                minInd = pt_j;
                        }
                }
                temp[pt_i] = minVal;
                done[minInd] = true;
        }
        for( int pt_i = 0; pt_i < N+1; pt_i++ ) { abcissa[pt_i] = temp[pt_i]; }

        for( int pt_i = 0; pt_i < N+1; pt_i++ ) {
                weights[pt_i] = 2.0/((1.0 - abcissa[pt_i]*abcissa[pt_i])*P->EvalDeriv(abcissa[pt_i])*P->EvalDeriv(abcissa[pt_i]));
        }

	for( int j = 0; j <= N; j++ ) {
		for( int i = 0; i <= N; i++ ) {
			xi[j*(N+1)+i] = new double[2];
			xi[j*(N+1)+i][0] = i*dx - 1.0;
			xi[j*(N+1)+i][1] = j*dx - 1.0;
			coord = new double[2];
			coord[0] = abcissa[i];
			coord[1] = abcissa[j];
			quadPts[j*(N+1)+i] = new QuadPoint( coord, weights[i]*weights[j] );
		}
	}

        delete[] temp;
        delete[] done;
	for( int poly_i = 0; poly_i < N+2; poly_i++ ) { delete polys[poly_i]; }
	delete[] polys;
}

Lagrange::~Lagrange() {}

#define EPS 1.0e-8

double Lagrange::Cj( int j, double xj, double x ) {
	double a, b;

	if( fabs( x - xj ) < EPS ) { return 1.0; }

	a = b = 1.0;
	for( int i = 0; i <= N; i++ ) {
		if( i == j ) { continue; }

		a *= ( x - xi[i][0] );
		b *= ( xj - xi[i][0] );
	}

	return a/b;
}

double* Lagrange::ShapeFuncs( int pt ) {
	return NULL;
}

double* Lagrange::ShapeFuncsAtCoord( double* coord ) {
	int i, j;

        for( int pt_i = 0; pt_i < nNodes; pt_i++ ) {
                i = pt_i%(N+1);
                j = pt_i/(N+1);
                Ni[pt_i] = Cj( i, xi[pt_i][0], coord[0] )*Cj( j, xi[pt_i][1], coord[1] );
        }

        return Ni;
}

double** Lagrange::ShapeFuncDerivs( int pt ) {
	//TODO
	return NULL;
}

double** Lagrange::ShapeFuncDerivsAtCoord( double* coord ) {
	//TODO
	return NULL;
}

/***********************************************************************************************************************/

/* note: assumes a single square element! */
Trig::Trig( int _N ) : Element( _N ) {
	int i, j;
	double* coord, dx;

	dx = 2.0/N;
	Ninv = 1.0/N;

	for( j = 0; j < N+1; j++ ) {
		for( i = 0; i < N+1; i++ ) {
			coord = new double[2];
			coord[0] = i*dx - 1.0;
			coord[1] = j*dx - 1.0;

			quadPts[j*(N+1)+i] = new QuadPoint( coord, 1 );
			xi[j*(N+1)+i] = new double[2];
			xi[j*(N+1)+i][0] = coord[0];
			xi[j*(N+1)+i][1] = coord[1];
		}
	}
}

Trig::~Trig() {}

double* Trig::ShapeFuncs( int pt ) {
	return NULL;
}

double* Trig::ShapeFuncsAtCoord( double* coord ) {
        int pt_i, i, j;

        for( pt_i = 0; pt_i < nNodes; pt_i++ ) {
                i = pt_i%(N+1);
                j = pt_i/(N+1);
		if( i == N || j == N ) {
			continue;
		}
                Ni[pt_i] = Cj( i, xi[pt_i][0], coord[0] )*Cj( j, xi[pt_i][1], coord[1] );
        }

	return Ni;
}

double** Trig::ShapeFuncDerivs( int pt ) {
	return NULL;
}

double** Trig::ShapeFuncDerivsAtCoord( double* coord ) {
	/* TODO */
	return NULL;
}

#define EPS 1.0e-8

double Trig::Cj( int j, double xj, double x ) {
	double a, b;

	if( fabs( x - xj ) < EPS ) {
		return 1.0;
	}
	if( j == 0 && fabs( x - 1.0 ) < EPS ) {
		return 1.0;
	}
	a = sin( 0.5*M_PI*N*(x - xj) );
	b = tan( 0.5*M_PI*(x - xj) );
	if( fabs( b ) < EPS ) {
		return 0.0;
	}
	return Ninv*a/b;
}

double Trig::dCjdx( int j, double xj, double x ) {
	/* TODO */
	return 0.0;
}

/***********************************************************************************************************************/

Linear::Linear( int _N ) : Element( _N ) {
	double* coord;

	weights[0] = 1.0;
	weights[1] = 1.0;
	abcissa[0] = -1.0/sqrt(3.0);
	abcissa[1] = +1.0/sqrt(3.0);

	coord = new double[2]; coord[0] = abcissa[0]; coord[1] = abcissa[0];
	quadPts[0] = new QuadPoint( coord, weights[0]*weights[0] );

	coord = new double[2]; coord[0] = abcissa[1]; coord[1] = abcissa[0];
	quadPts[1] = new QuadPoint( coord, weights[1]*weights[0] );

	coord = new double[2]; coord[0] = abcissa[0]; coord[1] = abcissa[1];
	quadPts[2] = new QuadPoint( coord, weights[0]*weights[1] );

	coord = new double[2]; coord[0] = abcissa[1]; coord[1] = abcissa[1];
	quadPts[3] = new QuadPoint( coord, weights[1]*weights[1] );

	xi[0] = new double[2]; xi[0][0] = -1.0;	xi[0][1] = -1.0;
	xi[1] = new double[2]; xi[1][0] = +1.0; xi[1][1] = -1.0;
	xi[2] = new double[2]; xi[2][0] = -1.0; xi[2][1] = +1.0;
	xi[3] = new double[2]; xi[3][0] = +1.0; xi[3][1] = +1.0;
}

Linear::~Linear() {}

double* Linear::ShapeFuncs( int pt ) {
	return ShapeFuncsAtCoord( quadPts[pt]->coord );
}

double* Linear::ShapeFuncsAtCoord( double* coord ) {
	double x = coord[0], y = coord[1];

        Ni[0] = 0.25*( 1.0-x )*( 1.0-y );
        Ni[1] = 0.25*( 1.0+x )*( 1.0-y );
        Ni[3] = 0.25*( 1.0+x )*( 1.0+y );
        Ni[2] = 0.25*( 1.0-x )*( 1.0+y );

	return Ni;
}

double** Linear::ShapeFuncDerivs( int pt ) {
	return ShapeFuncDerivsAtCoord( quadPts[pt]->coord );
}

double** Linear::ShapeFuncDerivsAtCoord( double* coord ) {
	double x = coord[0], y = coord[1];

        GNix[0][0] = - 0.25*( 1.0 - y );
        GNix[0][1] =   0.25*( 1.0 - y );
        GNix[0][3] =   0.25*( 1.0 + y );
        GNix[0][2] = - 0.25*( 1.0 + y );

        GNix[1][0] = - 0.25*( 1.0 - x );
        GNix[1][1] = - 0.25*( 1.0 + x );
        GNix[1][3] =   0.25*( 1.0 + x );
        GNix[1][2] =   0.25*( 1.0 - x );

	return GNix;
}

Quadratic::Quadratic( int _N ) : Element( _N ) {
	double* coord;

	weights[0] = 5.0/9.0;
	weights[1] = 8.0/9.0;
	weights[2] = 5.0/9.0;

	abcissa[0] = -sqrt(15.0)/5.0;
	abcissa[1] = 0.0;
	abcissa[2] = +sqrt(15.0)/5.0;

	coord = new double[2]; coord[0] = abcissa[0]; coord[1] = abcissa[0];
	quadPts[0] = new QuadPoint( coord, weights[0]*weights[0] );

	coord = new double[2]; coord[0] = abcissa[1]; coord[1] = abcissa[0];
	quadPts[1] = new QuadPoint( coord, weights[1]*weights[0] );

	coord = new double[2]; coord[0] = abcissa[2]; coord[1] = abcissa[0];
	quadPts[2] = new QuadPoint( coord, weights[2]*weights[0] );

	coord = new double[2]; coord[0] = abcissa[0]; coord[1] = abcissa[1];
	quadPts[3] = new QuadPoint( coord, weights[0]*weights[1] );

	coord = new double[2]; coord[0] = abcissa[1]; coord[1] = abcissa[1];
	quadPts[4] = new QuadPoint( coord, weights[1]*weights[1] );

	coord = new double[2]; coord[0] = abcissa[2]; coord[1] = abcissa[1];
	quadPts[5] = new QuadPoint( coord, weights[2]*weights[1] );

	coord = new double[2]; coord[0] = abcissa[0]; coord[1] = abcissa[2];
	quadPts[6] = new QuadPoint( coord, weights[0]*weights[2] );

	coord = new double[2]; coord[0] = abcissa[1]; coord[1] = abcissa[2];
	quadPts[7] = new QuadPoint( coord, weights[1]*weights[2] );

	coord = new double[2]; coord[0] = abcissa[2]; coord[1] = abcissa[2];
	quadPts[8] = new QuadPoint( coord, weights[2]*weights[2] );

	xi[0] = new double[2]; xi[0][0] = -1.0; xi[0][1] = -1.0;
	xi[1] = new double[2]; xi[1][0] =  0.0; xi[1][1] = -1.0;
	xi[2] = new double[2]; xi[2][0] = +1.0; xi[2][1] = -1.0;

	xi[3] = new double[2]; xi[3][0] = -1.0; xi[3][1] =  0.0;
	xi[4] = new double[2]; xi[4][0] =  0.0; xi[4][1] =  0.0;
	xi[5] = new double[2]; xi[5][0] = +1.0; xi[5][1] =  0.0;

	xi[6] = new double[2]; xi[6][0] = -1.0; xi[6][1] = +1.0;
	xi[7] = new double[2]; xi[7][0] =  0.0; xi[7][1] = +1.0;
	xi[8] = new double[2]; xi[8][0] = +1.0; xi[8][1] = +1.0;
}

Quadratic::~Quadratic() {}

double* Quadratic::ShapeFuncs( int pt ) {
	return ShapeFuncsAtCoord( quadPts[pt]->coord );
}

double* Quadratic::ShapeFuncsAtCoord( double* coord ) {
        double  x = coord[0], y = coord[1];
        double  a0 = x - 1.0, b0 = y - 1.0;
        double  a1 = 1.0 - x * x, b1 = 1.0 - y * y;
        double  a2 = x + 1.0, b2 = y + 1.0;
        double  m0 = 0.5 * x;
        double  m1 = 0.5 * y;
        double  m2 = 0.25 * x * y;

        Ni[0] = m2 * a0 * b0;
        Ni[1] = m1 * a1 * b0;
        Ni[2] = m2 * a2 * b0;

        Ni[3] = m0 * a0 * b1;
        Ni[4] = a1 * b1;
        Ni[5] = m0 * a2 * b1;

        Ni[6] = m2 * a0 * b2;
        Ni[7] = m1 * a1 * b2;
        Ni[8] = m2 * a2 * b2;

	return Ni;
}

double** Quadratic::ShapeFuncDerivs( int pt ) {
	return ShapeFuncDerivsAtCoord( quadPts[pt]->coord );
}

double** Quadratic::ShapeFuncDerivsAtCoord( double* coord ) {
        double  x = coord[0], y = coord[1];
        double  a0 = x - 1.0, b0 = y - 1.0;
        double  a1 = x + 1.0, b1 = y + 1.0;
        double  a2 = 2.0 * x - 1.0, b2 = 2.0 * y - 1.0;
        double  a3 = 2.0 * x + 1.0, b3 = 2.0 * y + 1.0;
        double  a4 = 1.0 - x * x, b4 = 1.0 - y * y;
        double  m0 = 0.25 * x;
        double  m1 = 0.25 * y;
        double  m2 = -x * y;

        /* Corner nodes. */
        GNix[0][0] = m1 * a2 * b0;
        GNix[0][2] = m1 * a3 * b0;
        GNix[0][6] = m1 * a2 * b1;
        GNix[0][8] = m1 * a3 * b1;
        GNix[1][0] = m0 * a0 * b2;
        GNix[1][2] = m0 * a1 * b2;
        GNix[1][6] = m0 * a0 * b3;
        GNix[1][8] = m0 * a1 * b3;

        /* Side nodes. */
        GNix[0][1] = m2 * b0;
        GNix[0][7] = m2 * b1;
        GNix[0][3] = 0.5 * a2 * b4;
        GNix[0][5] = 0.5 * a3 * b4;
        GNix[1][1] = 0.5 * a4 * b2;
        GNix[1][7] = 0.5 * a4 * b3;
        GNix[1][3] = m2 * a0;
        GNix[1][5] = m2 * a1;

        /* Center node. */
        GNix[0][4] = -2.0 * x * b4;
        GNix[1][4] = -2.0 * y * a4;

	return GNix;
}
