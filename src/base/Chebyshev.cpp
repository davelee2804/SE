#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Chebyshev.h"

using namespace std;

Chebyshev::Chebyshev( int _N ) : Element( _N ) {
	double* coord;

	for( int i = 0; i <= N; i++ ) {
		abcissa[i] = -cos( (M_PI*i)/N ); //in original definition abcissa go from -1=>+1
		weights[i] = M_PI/N;
	}
	weights[0] *= 0.5;
	weights[N] *= 0.5;

	for( int j = 0; j <= N; j++ ) {
		for( int i = 0; i <= N; i++ ) {
			coord = new double[2];
			coord[0] = abcissa[i];
			coord[1] = abcissa[j];
			quadPts[j*(N+1)+i] = new QuadPoint( coord, weights[i]*weights[j] );

			xi[j*(N+1)+i] = new double[2];
			xi[j*(N+1)+i][0] = abcissa[i];
			xi[j*(N+1)+i][1] = abcissa[j];
		}
	}

	/* for Chebychev polynomials, T_n = P_n^{alpha,beta}; alpha = beta = 0.5 */
	polys = new Jacobi*[N+1];
	for( int i = 0; i <= N; i++ ) {
		polys[i] = new Jacobi( i, 0.5 );
	}
	polys[0]->c[0] = 1.0;
	polys[1]->c[1] = 1.0;
	for( int i = 2; i <= N; i++ ) {
		for( int j = 0; j <= i-2; j++ ) {
			polys[i]->c[j] -= polys[i-2]->c[j];
		}
		for( int j = 1; j <= i; j++ ) {
			polys[i]->c[j] += 2*polys[i-1]->c[j-1];
		}
	}

	P = new Jacobi( N, 0.5 );
	for( int j = 0; j <= N-2; j++ ) {
		P->c[j] -= polys[N-2]->c[j];
	}
	for( int j = 1; j <= N; j++ ) {
		P->c[j] += 2*polys[N-1]->c[j-1];
	}

	dCij     = new double[(N+1)*(N+1)];
	d2Cij    = new double[(N+1)*(N+1)];
	GNixx    = new double*[3];
	GNixx[0] = new double[nNodes];
	GNixx[1] = new double[nNodes];
	GNixx[2] = new double[nNodes];

	for( int i = 0; i <= N; i++ ) {
		for( int j = 0; j <= N; j++ ) {
			/* switch these two as abcissa go from -1 to +1 */
			if( i == 0 && j == 0 ) {
				dCij[i*(N+1)+j] = -(1.0 + 2.0*N*N)/6.0;
			}
			else if( i == N && j == N ) {
				dCij[i*(N+1)+j] = +(1.0 + 2.0*N*N)/6.0;
			}
			else if( i == j ) {
				dCij[i*(N+1)+j] = -0.5*abcissa[i]/(1.0 - abcissa[i]*abcissa[i]);
			}
			else {
				double pi = (i == 0 || i == N) ? 2.0 : 1.0;
				double pj = (j == 0 || j == N) ? 2.0 : 1.0;
				dCij[i*(N+1)+j] = pow(-1,i+j)*pi/(pj*(abcissa[i] - abcissa[j]));
			}
			//dCij[i*(N+1)+j] = dCjdx(j,abcissa[i]);
		}
	}
	Mult( dCij, dCij, d2Cij, N+1 );

/*
*/
for(int i=0;i<=N;i++){
cout<<P->c[i]<<endl;
}
char filename[40];
ofstream file;
double x;
for(int i=0;i<=N;i++){
  sprintf(filename,"c_%.2u.txt",i);
  file.open(filename);
  for(int j=0;j<=200;j++){
    x=-cos(M_PI*j/200.0);
    file<<Cj(i,x)<<endl;
  }
  file.close();
}
for(int i=0;i<=N;i++){
  sprintf(filename,"dc_%.2u.txt",i);
  file.open(filename);
  file<<dCjdx(i,-1.0)<<"\t"<<(Cj(i,-1.0+0.01)-Cj(i,-1.0))/0.01<<endl;
  for(int j=1;j<=199;j++){
    x=-1.0+j*0.01;
    file<<dCjdx(i,x)<<"\t"<<(Cj(i,x+0.01)-Cj(i,x-0.01))/0.02<<endl;
  }
  file<<dCjdx(i,+1.0)<<"\t"<<(Cj(i,+1.0)-Cj(i,+1.0-0.01))/0.01<<endl;
  file.close();
}
for(int i=0;i<=N;i++){
  sprintf(filename,"dCij_%.2u.txt",i);
  file.open(filename);
  for(int j=0;j<=N;j++){
    //file<<dCij[i*(N+1)+j]<<endl;
    file<<dCij[j*(N+1)+i]<<endl;
  }
  file.close();
}
}

Chebyshev::~Chebyshev() {
	delete[] dCij;
	delete[] d2Cij;
	delete[] GNixx[0];
	delete[] GNixx[1];
	delete[] GNixx[2];
	delete[] GNixx;
	for( int i = 0; i <= N; i++ ) {
		delete polys[i];
	}
	delete[] polys;
}

double* Chebyshev::ShapeFuncs( int pt ) {
	for( int i = 0; i < nNodes; i++ ) {
		Ni[i] = 0.0;
	}
	Ni[pt] = 1.0;

	return Ni;
}

double* Chebyshev::ShapeFuncsAtCoord( double* x ) {
        int pt_i, i, j;

        for( pt_i = 0; pt_i < nNodes; pt_i++ ) {
                i = pt_i%(N+1);
                j = pt_i/(N+1);
                Ni[pt_i] = Cj( i, x[0] )*Cj( j, x[1] );
        }

        return Ni;
}

double** Chebyshev::ShapeFuncDerivs( int pt ) {
	int x_i, y_i, x_j, y_j;

	x_i = pt%(N+1);
	y_i = pt/(N+1);

	for( int pt_j = 0; pt_j < nNodes; pt_j++ ) {
		x_j = pt_j%(N+1);
		y_j = pt_j/(N+1);
		GNix[0][pt_j] = (y_i == y_j) ? dCij[x_i*(N+1)+x_j] : 0.0;
		GNix[1][pt_j] = (x_i == x_j) ? dCij[y_i*(N+1)+y_j] : 0.0;
	}
	return GNix;
}

double** Chebyshev::ShapeFuncDerivsAtCoord( double* x ) {
        int pt_i, i, j;

        for( pt_i = 0; pt_i < nNodes; pt_i++ ) {
                i = pt_i%(N+1);
                j = pt_i/(N+1);
                GNix[0][pt_i] = dCjdx( i, x[0] )*Cj( j, x[1] );
                GNix[1][pt_i] = Cj( i, x[0] )*dCjdx( j, x[1] );
        }

        return GNix;
}

double** Chebyshev::ShapeFunc2ndDerivs( int pt ) {
	int x_i, y_i, x_j, y_j;

	x_i = pt%(N+1);
	y_i = pt/(N+1);

	for( int pt_j = 0; pt_j < nNodes; pt_j++ ) {
		x_j = pt_j%(N+1);
		y_j = pt_j/(N+1);
		GNixx[0][pt_j] = (y_i == y_j) ? d2Cij[x_i*(N+1)+x_j] : 0.0;
		GNixx[1][pt_j] = (x_i == x_j) ? d2Cij[y_i*(N+1)+y_j] : 0.0;
		GNixx[2][pt_j] = dCij[x_i*(N+1)+x_j]*dCij[y_i*(N+1)+y_j];
	}
	return GNixx;
}

double Chebyshev::Cj( int j, double x ) {
        double cj = 0.0, cm, xj = abcissa[j];

	for( int m = 0; m <= N; m++ ) {
		cm = (polys[m]->Eval(x))*(polys[m]->Eval(xj));
		if( m == 0 || m == N ) {
			cm *= 0.5;
		}
		cj += cm;
	}
	cj *= 2.0/N;
	if( j == 0 || j == N ) {
		cj *= 0.5;
	}

        return cj;
}

double Chebyshev::dCjdx( int j, double x ) {
        double dcj = 0.0, dcm, xj = abcissa[j];

	for( int m = 0; m <= N; m++ ) {
		dcm = (polys[m]->EvalDeriv(x))*(polys[m]->Eval(xj));
		if( m == 0 || m == N ) {
			dcm *= 0.5;
		}
		dcj += dcm;
	}
	dcj *= 2.0/N;
	if( j == 0 || j == N ) {
		dcj *= 0.5;
	}

        return dcj;
}
