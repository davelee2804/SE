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
#include "NeumannBC.h"

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

NeumannBC::NeumannBC( string _name, Mesh* _mesh, Field* _velocity, Field* _prevVel, Field* _velHat, double _nu ) : RHSOp( _name, _mesh, 0.0, NULL ) {
	Field* 		currVel;
	int 		node_i;
	double**	dv;
        int 		pt_i, pt_j, endLoop;
        double 		xCurr, xPrev, x0;
        bool* 		done;
        double* 	temp;
        double 		minVal;
        int 		minInd;
        Jacobi* 	P 		= mesh->el->polys[mesh->el->N+1];

	velocity = _velocity;
	prevVel  = _prevVel;
	velHat   = _velHat;
	nu       = _nu;

	/* set up the (extrapolated) current time step vorticity field */
	currVel = new Field( "currVel", velocity->mesh, velocity->nDofs, NULL );
	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		currVel->vals[node_i][0] = 2.0*velocity->vals[node_i][0] - prevVel->vals[node_i][0];
		currVel->vals[node_i][1] = 2.0*velocity->vals[node_i][1] - prevVel->vals[node_i][1];
	}

	vorticity = new Field( "vorticity", velocity->mesh, 1, NULL );

	dv = new double*[2];
	dv[0] = new double[2];
	dv[1] = new double[2];
	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		currVel->InterpDerivsGlobal( velocity->mesh->verts[node_i], dv );
		vorticity->vals[node_i][0] = dv[1][0] - dv[0][1];
	}
	delete[] dv[0];
	delete[] dv[1];
	delete dv;

	delete currVel;

	/* calculate the gaussian quadrature weights and abcissa */
        weight  = new double[mesh->el->N+1];
        abcissa = new double[mesh->el->N+1];
        endLoop = (mesh->el->N+1)%2 == 1 ? mesh->el->N/2 : (mesh->el->N+1)/2;
        for( pt_i = 0; pt_i < endLoop; pt_i++ ) {
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
        for( pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
                done[pt_i] = false;
        }
        for( pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
                minVal = 1.0e+99;
                for( pt_j = 0; pt_j < mesh->el->N+1; pt_j++ ) {
                        if( !done[pt_j] && abcissa[pt_j] < minVal ) {
                                minVal = abcissa[pt_j];
                                minInd = pt_j;
                        }
                }
                temp[pt_i] = minVal;
                done[minInd] = true;
        }
        for( pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
                abcissa[pt_i] = temp[pt_i];
        }
        delete[] temp;
        delete[] done;

        for( pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
                weight[pt_i] = 2.0/((1.0 - abcissa[pt_i]*abcissa[pt_i])*P->EvalDeriv(abcissa[pt_i])*P->EvalDeriv(abcissa[pt_i]));
        }

        dwdx = new double*[1];
        dwdx[0] = new double[2];
        Ni = new double[mesh->el->N+1];
}

NeumannBC::~NeumannBC() {
        delete[] weight;
        delete[] abcissa;
        delete[] dwdx[0];
        delete[] dwdx;
        delete[] Ni;
	delete vorticity;
}

void NeumannBC::AssembleElement( int el_i, double* rhs ) {
	Element* 	el = mesh->el;
	double 		detJac, gCoordBL[2], gCoordTR[2], gCoord[2], xi[2], v[2], norm, a;
	int		point_i, node_i, offset, stride;
	double*		Ni;

	mesh->LocalToGlobal( el->quadPts[0]->coord, el_i, gCoordBL );
	mesh->LocalToGlobal( el->quadPts[el->nPoints-1]->coord, el_i, gCoordTR );

	/* bottom */
	if( gCoordBL[1] < mesh->min[1] + EPS ) {
		detJac = (mesh->max[0] - mesh->min[0])/(2.0*mesh->nEls[0]);
		norm = -1.0;
		for( point_i = 0; point_i < el->N+1; point_i++ ) {
			a = detJac*weight[point_i]*nu*norm;
			Ni = ShapeFuncs( abcissa[point_i] );
			xi[0] = abcissa[point_i];
			xi[1] = -1.0;
			velHat->InterpLocal( el_i, xi, v );
			mesh->LocalToGlobal( xi, el_i, gCoord );
			vorticity->InterpDerivsGlobal( gCoord, dwdx );
			for( node_i = 0; node_i < el->N+1; node_i++ ) {
				rhs[node_i] -= a*( v[1] - dwdx[0][0] );
			}
		}
	}
	/* top */
	if( gCoordTR[1] > mesh->max[1] - EPS ) {
		detJac = (mesh->max[0] - mesh->min[0])/(2.0*mesh->nEls[0]);
		norm = +1.0;
		offset = (el->N+1)*el->N;
		for( point_i = 0; point_i < el->N+1; point_i++ ) {
			a = detJac*weight[point_i]*nu*norm;
			Ni = ShapeFuncs( abcissa[point_i] );
			xi[0] = abcissa[point_i];
			xi[1] = +1.0;
			velHat->InterpLocal( el_i, xi, v );
			mesh->LocalToGlobal( xi, el_i, gCoord );
			vorticity->InterpDerivsGlobal( gCoord, dwdx );
			for( node_i = 0; node_i < el->N+1; node_i++ ) {
				rhs[offset+node_i] -= a*( v[1] - dwdx[0][0] );
			}
		}
	}
	/* left */
	if( gCoordBL[0] < mesh->min[0] + EPS ) {
		detJac = (mesh->max[1] - mesh->min[1])/(2.0*mesh->nEls[1]);
		norm = -1.0;
		stride = el->N+1;
		for( point_i = 0; point_i < el->N+1; point_i++ ) {
			a = detJac*weight[point_i]*nu*norm;
			Ni = ShapeFuncs( abcissa[point_i] );
			xi[0] = -1.0;
			xi[1] = abcissa[point_i];
			velHat->InterpLocal( el_i, xi, v );
			mesh->LocalToGlobal( xi, el_i, gCoord );
			vorticity->InterpDerivsGlobal( gCoord, dwdx );
			for( node_i = 0; node_i < el->N+1; node_i++ ) {
				rhs[node_i*stride] -= a*( v[0] + dwdx[0][1] );
			}
		}
	}
	/* right */
	if( gCoordTR[0] > mesh->max[0] - EPS ) {
		detJac = (mesh->max[1] - mesh->min[1])/(2.0*mesh->nEls[1]);
		norm = +1.0;
		stride = el->N+1;
		offset = el->N;
		for( point_i = 0; point_i < el->N+1; point_i++ ) {
			a = detJac*weight[point_i]*nu*norm;
			Ni = ShapeFuncs( abcissa[point_i] );
			xi[0] = +1.0;
			xi[1] = abcissa[point_i];
			velHat->InterpLocal( el_i, xi, v );
			mesh->LocalToGlobal( xi, el_i, gCoord );
			vorticity->InterpDerivsGlobal( gCoord, dwdx );
			for( node_i = 0; node_i < el->N+1; node_i++ ) {
				rhs[offset+(node_i*stride)] -= a*( v[0] + dwdx[0][1] );
			}
		}
	}
}

double* NeumannBC::ShapeFuncs( double x ) {
        //int i;
	//Legendre* el = (Legendre*)mesh->el;

/*
        for( i = 0; i <= el->N; i++ ) {
                Ni[i] = el->Cj( i, mesh->el->xi[i][0], x );
        }
        return Ni;
*/
	return NULL;
}
