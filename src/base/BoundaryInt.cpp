#include <iostream>
#include <string>
#include <cmath>

#include "QuadPoint.h"
#include "Jacobi.h"
#include "Element.h"
#include "Legendre.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "RHSOp.h"
#include "BoundaryInt.h"

using namespace std;
using std::string;

#define STEP 1.0e-3
#define TOL 1.0e-12
#define EPS 1.0e-6

bool Found( int n, double* roots, double root ) {
        int i;

        for( i = 0; i < n; i++ ) {
                if( fabs(roots[i] - root) < 1.0e-6 ) {
                        return true;
                }
        }
        return false;
}

BoundaryInt::BoundaryInt( string _name, Mesh* _mesh, Field* _field, Field* _fPrev ) : RHSOp( _name, _mesh, 0.0, _field ) {
        int 		endLoop;
        double 		xCurr, xPrev, x0;
        bool* 		done;
        double* 	temp;
        double 		minVal;
        int 		minInd;
        Jacobi* 	P 		= mesh->el->polys[mesh->el->N+1];

	fPrev = _fPrev;

	/* calculate the gaussian quadrature weights and abcissa */
        weight  = new double[mesh->el->N+1];
        abcissa = new double[mesh->el->N+1];
        endLoop = (mesh->el->N+1)%2 == 1 ? mesh->el->N/2 : (mesh->el->N+1)/2;
        for( int pt_i = 0; pt_i < endLoop; pt_i++ ) {
                x0 = (pt_i == 0) ? -1.0 : abcissa[pt_i-1];
                do {
                        x0 += STEP;
                        xCurr = x0;
                        do {
                                xPrev = xCurr;
                                xCurr = xPrev - P->Eval( xPrev )/P->EvalDeriv( xPrev );
                        } while( fabs( xCurr - xPrev ) > TOL );
                } while( Found( pt_i, abcissa, xCurr ) || xCurr > 0.0 - STEP );
                abcissa[pt_i] = xCurr;
                abcissa[mesh->el->N-pt_i] = -abcissa[pt_i];
        }
        if( (mesh->el->N+1)%2 == 1 ) {
                abcissa[(mesh->el->N+1)/2] = 0.0;
        }
        done = new bool[mesh->el->N+1];
        temp = new double[mesh->el->N+1];
        minInd = 999999999;
        for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
                done[pt_i] = false;
        }
        for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
                minVal = 1.0e+99;
                for( int pt_j = 0; pt_j < mesh->el->N+1; pt_j++ ) {
                        if( !done[pt_j] && abcissa[pt_j] < minVal ) {
                                minVal = abcissa[pt_j];
                                minInd = pt_j;
                        }
                }
                temp[pt_i] = minVal;
                done[minInd] = true;
        }
        for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
                abcissa[pt_i] = temp[pt_i];
        }
        delete[] temp;
        delete[] done;

        for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
                weight[pt_i] = 2.0/((1.0 - abcissa[pt_i]*abcissa[pt_i])*P->EvalDeriv(abcissa[pt_i])*P->EvalDeriv(abcissa[pt_i]));
        }
}

BoundaryInt::~BoundaryInt() {
        delete[] weight;
        delete[] abcissa;
}

/* Assumes use of Legendre type shape funcs, which are orthogonal to the quadrature 
   points so only evaluate the value at the specific quadrature point. */

void BoundaryInt::AssembleElement( int el_i, double* rhs ) {
	Element* 	el = mesh->el;
	double 		detJac, gCoordBL[2], gCoordTR[2], xi[2], norm, a, vCurr, vPrev;
	int		point_i, node_i, offset, stride;

	mesh->LocalToGlobal( el->quadPts[0]->coord, el_i, gCoordBL );
	mesh->LocalToGlobal( el->quadPts[el->nPoints-1]->coord, el_i, gCoordTR );

	/* bottom */
	if( !field->bcs->bcBottom[0] && !field->mesh->periodic[1] ) {
		if( gCoordBL[1] < mesh->min[1] + EPS ) {
			detJac = (mesh->max[0] - mesh->min[0])/(2.0*mesh->nEls[0]);
			norm = -1.0;
			for( point_i = 0; point_i < el->N+1; point_i++ ) {
				a = detJac*weight[point_i]*norm;
				xi[0] = abcissa[point_i];
				xi[1] = -1.0;
				field->InterpLocal( el_i, xi, &vCurr );
				fPrev->InterpLocal( el_i, xi, &vPrev );
				node_i = point_i;
				rhs[node_i] -= a*(2*vCurr - vPrev);
			}
		}
	}
	/* top */
	if( !field->bcs->bcTop[0] && !field->mesh->periodic[1] ) {
		if( gCoordTR[1] > mesh->max[1] - EPS ) {
			detJac = (mesh->max[0] - mesh->min[0])/(2.0*mesh->nEls[0]);
			norm = +1.0;
			offset = (el->N+1)*el->N;
			for( point_i = 0; point_i < el->N+1; point_i++ ) {
				a = detJac*weight[point_i]*norm;
				xi[0] = abcissa[point_i];
				xi[1] = +1.0;
				field->InterpLocal( el_i, xi, &vCurr );
				fPrev->InterpLocal( el_i, xi, &vPrev );
				node_i = point_i;
				rhs[offset+node_i] -= a*(2*vCurr - vPrev);
			}
		}
	}
	/* left */
	if( !field->bcs->bcLeft[0] && !field->mesh->periodic[0] ) {
		if( gCoordBL[0] < mesh->min[0] + EPS ) {
			detJac = (mesh->max[1] - mesh->min[1])/(2.0*mesh->nEls[1]);
			norm = -1.0;
			stride = el->N+1;
			for( point_i = 0; point_i < el->N+1; point_i++ ) {
				a = detJac*weight[point_i]*norm;
				xi[0] = -1.0;
				xi[1] = abcissa[point_i];
				field->InterpLocal( el_i, xi, &vCurr );
				fPrev->InterpLocal( el_i, xi, &vPrev );
				node_i = point_i;
				rhs[node_i*stride] -= a*(2*vCurr - vPrev);
			}
		}
	}
	/* right */
	if( !field->bcs->bcRight[0] && !field->mesh->periodic[0] ) {
		if( gCoordTR[0] > mesh->max[0] - EPS ) {
			detJac = (mesh->max[1] - mesh->min[1])/(2.0*mesh->nEls[1]);
			norm = +1.0;
			stride = el->N+1;
			offset = el->N;
			for( point_i = 0; point_i < el->N+1; point_i++ ) {
				a = detJac*weight[point_i]*norm;
				xi[0] = +1.0;
				xi[1] = abcissa[point_i];
				field->InterpLocal( el_i, xi, &vCurr );
				fPrev->InterpLocal( el_i, xi, &vPrev );
				node_i = point_i;
				rhs[offset+(node_i*stride)] -= a*(2*vCurr - vPrev);
			}
		}
	}
}
