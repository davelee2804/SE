#include <cmath>
#include <iostream>

#include "LinAlg.h"

#define I(p,q,N) \
        p*N + q

using namespace std;
using std::string;

/* A^T_{j,i} = A_{i,j} */
void Tran( double* A, double* AT, int n ) {
        int i, j;
        for( i = 0; i < n; i++ ) {
                for( j = 0; j < n; j++ ) {
                        AT[I(j,i,n)] = A[I(i,j,n)];
                }
        }
}

/* C_{nm} = A_{nl}*B_{lm} */
void Mult( double* A, double* B, double* C, int n ) {
        int i, j, k;
        for( i = 0; i < n*n; i++ ) {
                C[i] = 0.0;
        }
        for( i = 0; i < n; i++ ) {
                for( j = 0; j < n; j++ ) {
                        for( k = 0; k < n; k++ ) {
                                C[I(i,j,n)] += A[I(i,k,n)]*B[I(k,j,n)];
                        }
                }
        }
}

/* b_i = A_{ij}x_j */
void AXEB( double* A, double* x, double* b, int n ) {
        int i, j;

        for( i = 0; i < n; i++ ) {
                b[i] = 0.0;
                for( j = 0; j < n; j++ ) {
                        b[i] += A[I(i,j,n)]*x[j];
                }
        }

}

/* Gauss-Jordan matrix inverse routine - see Numerical Recipes in C, 2nd Ed.section 2.1 for details */
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void Inv( double* A, double* Ainv, int n ) {
        int *indxc, *indxr, *ipiv;
        int i, j, k, l, irow = 0, icol = 0, ll;
        double big, dum, pivinv, temp;

        indxc = new int[n]; indxr = new int[n]; ipiv  = new int[n];

        for( i = 0; i < n*n; i++ ) { Ainv[i] = A[i]; }
        for( j = 0; j< n; j++ ) { ipiv[j] = 0; }
        for( i = 0; i < n; i++ ) {
                big = 0.0;
                for( j = 0; j < n; j++ ) {
                        if( ipiv[j] != 1 ) {
                                for( k = 0; k < n; k++ ) {
                                        if( ipiv[k] == 0 ) {
                                                if( fabs(Ainv[I(j,k,n)]) >= big ) {
                                                        big = fabs(Ainv[I(j,k,n)]);
                                                        irow = j;
                                                        icol = k;
                                                }
                                        }
                                        else if( ipiv[k] > 1 ) { cerr << "Matrix inverse error! - singular matrix (1)\n"; }
                                }
                        }
                }
                ++(ipiv[icol]);
                if( irow != icol ) {
                        for( l = 0; l < n; l++ ) {
                                SWAP( Ainv[I(irow,l,n)], Ainv[I(icol,l,n)] );
                        }
                }
                indxr[i] = irow;
                indxc[i] = icol;
                if( fabs(Ainv[I(icol,icol,n)]) < 1.0e-12 ) { cerr << "Matrix inverse error! - singular matrix (2)\n"; }
                pivinv = 1.0/Ainv[I(icol,icol,n)];
                Ainv[I(icol,icol,n)] = 1.0;
                for( l = 0; l < n; l++ ) { Ainv[I(icol,l,n)] *= pivinv; }
                for( ll = 0; ll < n; ll++ ) {
                        if( ll != icol ) {
                                dum = Ainv[I(ll,icol,n)];
                                Ainv[I(ll,icol,n)] = 0.0;
                                for( l = 0; l < n; l++ ) { Ainv[I(ll,l,n)] -= Ainv[I(icol,l,n)]*dum; }
                        }
                }
        }
        for( l = n-1; l >= 0; l-- ) {
                if( indxr[l] != indxc[l] ) {
                        for( k = 0; k < n; k++ ) {
                                SWAP( Ainv[I(k,indxr[l],n)], Ainv[I(k,indxc[l],n)] );
                        }
                }
        }
        delete[] indxc; delete[] indxr; delete[] ipiv;
}
