#include <string>
#include <iostream>
#include <cstdio>
#include <cmath>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Advector.h"

using namespace std;
using std::string;

/* Second order accurate semi-Lagrangian advection scheme for spectral elements. 
 * References: 
 * 	Giraldo, F. X., J. B. Perot and P. F. Fischer (2003), "A spectral element semi-Lagrangian 
 * 	(SESL) method for the spherical shallow water equations", Journal of Computational 
 * 	Physics, 190, 623-650 
 *	Xiu, D & G. E. Karniadakis (2001) "A Semi-Lagrangian High-Order Method for Navier-Stokes
 *	Equations" Journal of Computational Physics, 172, 658-684
 */

Advector::Advector( Field* _field, Field* _velocity ) {
	char fieldName[40];

	field = _field;
	velocity = _velocity;

	secondOrder = false;
	cubicInterp = (field->mesh->N < 3) ? true : false;

	sprintf( fieldName, "%s-SL", field->name.c_str() );
	fieldSL = new Field( fieldName, field->mesh, field->nDofs, field->bcs );
	fieldSL->Copy( field );

	ptsX = new double*[4];
	ptsX[0] = new double[field->nDofs];
	ptsX[1] = new double[field->nDofs];
	ptsX[2] = new double[field->nDofs];
	ptsX[3] = new double[field->nDofs];
	ptsY = new double*[4];
	ptsY[0] = new double[field->nDofs];
	ptsY[1] = new double[field->nDofs];
	ptsY[2] = new double[field->nDofs];
	ptsY[3] = new double[field->nDofs];
}

Advector::Advector( Field* _field, Field* _velocity, Field* _fieldMinusOne, Field* _velMinusOne ) {
	char fieldName[40];

	field = _field;
	velocity = _velocity;
	fieldMinusOne = _fieldMinusOne;
	velMinusOne = _velMinusOne;

	secondOrder = true;
	cubicInterp = (field->mesh->N < 3) ? true : false;

	sprintf( fieldName, "%s-SL", field->name.c_str() );
	fieldSL = new Field( fieldName, field->mesh, field->nDofs, field->bcs );	
	fieldSL->Copy( field ); /* copy the bcs over to the semi-Lagrangian fields */

	sprintf( fieldName, "%s_n+1", velocity->name.c_str() );
	velPlusOne = new Field( fieldName, velocity->mesh, velocity->nDofs, velocity->bcs );

	sprintf( fieldName, "%s_n+1/2", velocity->name.c_str() );
	velPlusHalf = new Field( fieldName, velocity->mesh, velocity->nDofs, velocity->bcs );

	sprintf( fieldName, "%s_n-1/2", velocity->name.c_str() );
	velMinusHalf = new Field( fieldName, velocity->mesh, velocity->nDofs, velocity->bcs );

	ptsX = new double*[4];
	ptsX[0] = new double[field->nDofs];
	ptsX[1] = new double[field->nDofs];
	ptsX[2] = new double[field->nDofs];
	ptsX[3] = new double[field->nDofs];
	ptsY = new double*[4];
	ptsY[0] = new double[field->nDofs];
	ptsY[1] = new double[field->nDofs];
	ptsY[2] = new double[field->nDofs];
	ptsY[3] = new double[field->nDofs];
}

Advector::~Advector() {
	delete fieldSL;

	if( secondOrder ) {
		delete velPlusOne;
		delete velPlusHalf;
		delete velMinusHalf;
	}

	delete[] ptsX[0];
	delete[] ptsX[1];
	delete[] ptsX[2];
	delete[] ptsX[3];
	delete[] ptsX;
	delete[] ptsY[0];
	delete[] ptsY[1];
	delete[] ptsY[2];
	delete[] ptsY[3];
	delete[] ptsY;
}

void Advector::Advect( double dt ) {
	cout << "advecting " << field->name.c_str() << " field...";

	if( secondOrder ) {
		ExtrapolateVelocities();
		SecondOrder( dt );
	}
	else {
		FirstOrder( dt );
	}

	cout << " done.\n";
}

void Advector::PeriodicUpdate( double* coord ) {
        Mesh* mesh = field->mesh;
        int dim_i;

        for( dim_i = 0; dim_i < field->mesh->dim; dim_i++ ) {
                if( coord[dim_i] < mesh->min[dim_i] ) {
                        coord[dim_i] = (mesh->periodic[dim_i]) ? mesh->max[dim_i] - mesh->min[dim_i] + coord[dim_i] : mesh->min[dim_i];
                }
                if( coord[dim_i] > mesh->max[dim_i] ) {
                        coord[dim_i] = (mesh->periodic[dim_i]) ? mesh->min[dim_i] - mesh->max[dim_i] + coord[dim_i] : mesh->max[dim_i];
                }
        }
}

void Advector::IntegrateRK2( double dt, Field* velField, double* origin, double* coord ) {
	double xHat[3], vel[3];

	velocity->InterpGlobal( origin, vel );
	xHat[0] = origin[0] - 0.5*dt*vel[0];
	xHat[1] = origin[1] - 0.5*dt*vel[1];
	PeriodicUpdate( xHat );

	velField->InterpGlobal( xHat, vel );
	coord[0] = origin[0] - dt*vel[0];
	coord[1] = origin[1] - dt*vel[1];
	PeriodicUpdate( coord );
}

void Advector::IntegrateRK4( double dt, double* origin, double* coord ) {
        double k[4][3];
        double coordPrime[3];
        int dim_i;

        velocity->InterpGlobal( origin, k[0] );

        for( dim_i = 0; dim_i < field->mesh->dim; dim_i++ ) {
                coordPrime[dim_i] = origin[dim_i] - 0.5*dt*k[0][dim_i];
        }
        PeriodicUpdate( coordPrime );
        velocity->InterpGlobal( coordPrime, k[1] );

        for( dim_i = 0; dim_i < field->mesh->dim; dim_i++ ) {
                coordPrime[dim_i] = origin[dim_i] - 0.5*dt*k[1][dim_i];
        }
        PeriodicUpdate( coordPrime );
        velocity->InterpGlobal( coordPrime, k[2] );

        for( dim_i = 0; dim_i < field->mesh->dim; dim_i++ ) {
                coordPrime[dim_i] = origin[dim_i] - dt*k[2][dim_i];
        }
        PeriodicUpdate( coordPrime );
        velocity->InterpGlobal( coordPrime, k[3] );

        for( dim_i = 0; dim_i < field->mesh->dim; dim_i++ ) {
                coord[dim_i] = origin[dim_i] - dt*(1.0/6.0)*(k[0][dim_i] + 2.0*k[1][dim_i] + 2.0*k[2][dim_i] + k[3][dim_i]);
        }
        PeriodicUpdate( coord );
}

void Advector::FirstOrder( double dt ) {
	double coord[3], val[3];
	int node_i, dof_i;

	for( node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
		IntegrateRK2( dt, velocity, field->mesh->verts[node_i], coord );
		if( cubicInterp ) {
			Interpolate( field, coord, val );
		}
		else {
			field->InterpGlobal( coord, val );
		}
		for( dof_i = 0; dof_i < field->nDofs; dof_i++ ) {
			fieldSL->vals[node_i][dof_i] = val[dof_i];
		}
	}

	if( field->mesh->periodic[0] || field->mesh->periodic[1] ) {
		fieldSL->PeriodicUpdate();
	}
}

void Advector::SecondOrder( double dt ) {
	int node_i, dof_i;
	double *origin, x1[3], x2[3], psi1[3], psi2[3], v[2], xHalf[2];

	for( node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
		origin = field->mesh->verts[node_i];

		velPlusOne->InterpGlobal( origin, v );
		xHalf[0] = origin[0] - 0.5*dt*v[0];
		xHalf[1] = origin[1] - 0.5*dt*v[1];
		PeriodicUpdate( xHalf );

		velPlusHalf->InterpGlobal( xHalf, v );
		x1[0] = origin[0] - dt*v[0];
		x1[1] = origin[1] - dt*v[1];
		PeriodicUpdate( x1 );

		velocity->InterpGlobal( x1, v );
		xHalf[0] = x1[0] - 0.5*dt*v[0];
		xHalf[1] = x1[1] - 0.5*dt*v[1];
		PeriodicUpdate( xHalf );

		velMinusHalf->InterpGlobal( xHalf, v );
		x2[0] = x1[0] - dt*v[0];
		x2[1] = x1[1] - dt*v[1];
		PeriodicUpdate( x2 );

		if( cubicInterp ) {
			Interpolate( field, x1, psi1 );
			Interpolate( fieldMinusOne, x2, psi2 );
		}
		else {
			field->InterpGlobal( x1, psi1 );
			fieldMinusOne->InterpGlobal( x2, psi2 );
		}
		for( dof_i = 0; dof_i < field->nDofs; dof_i++ ) {
			fieldSL->vals[node_i][dof_i] = 2.0*psi1[dof_i] - 0.5*psi2[dof_i];
		}
	}

	if( field->mesh->periodic[0] || field->mesh->periodic[1] ) {
		fieldSL->PeriodicUpdate();
	}
}

void Advector::InterpLagrange( double x, double* coords, int nDofs, double** vals, double* psi ) {
	int 	node_i, dof_i;
	int	otherIndices[3];
        int	otherIndexCount, otherIndex_i;
        double	factor;

	for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
		psi[dof_i] = 0.0;
	}

	for( node_i = 0; node_i < 4; node_i++ ) {
		otherIndexCount = 0;
		for( otherIndex_i = 0; otherIndex_i < 4; otherIndex_i++ ) {
			if( otherIndex_i != node_i ) {
				otherIndices[otherIndexCount++] = otherIndex_i;
			}
		}

                factor = 1.0;
		for( otherIndex_i = 0; otherIndex_i < 3; otherIndex_i++ ) {
			factor *= (x - coords[otherIndices[otherIndex_i]])/(coords[node_i] - coords[otherIndices[otherIndex_i]]);
		}

		for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
			psi[dof_i] += vals[node_i][dof_i]*factor;
		}
        }
}

/* bicubic interpolation scheme - assumes that nodes are evenly spaced */
void Advector::Interpolate( Field* field, double* coord, double* psi ) {
	Mesh*		mesh		= field->mesh;	
	int 		*nodes, el_i, node, dof_i;
	int		x_0, y_0;
	int		x_i, y_i;
	int		nodeIndex[4][4];
	double		px[4], py[4];
	double		lCoord[2];
	double		dx = mesh->dx[0]/mesh->N;
	double		dy = mesh->dx[1]/mesh->N;

	mesh->GlobalToLocal( coord, &el_i, lCoord );
	nodes = mesh->ElNodes( el_i );
	node = nodes[0];

	x_0 = (int) node % mesh->nVerts[0];
	y_0 = ( (int) node / mesh->nVerts[0] ) % mesh->nVerts[1];

	/* get the number of nodes across and up that the point lies... */
	/* bottom left corner of stencil is closer to LHS of element */
	if( coord[0] <= mesh->verts[nodes[1]][0] ) x_0--;
	/* bottom left corner of stencil is closer to bottom of element */
	if( coord[1] <= mesh->verts[nodes[3]][1] ) y_0--;
	
	/* LHS node is global domain boundary */
	if( coord[0] <= mesh->min[0] + dx ) x_0++;
	/* RHS node is global domain boundary */
	else if( coord[0] >= mesh->max[0] - dx ) x_0--;

	/* top node is global domain boundary */
	if( coord[1] <= mesh->min[1] + dy ) y_0++;
	/* bottom node is global domain boundary */
	else if( coord[1] >= mesh->max[1] - dy ) y_0--;

	/* interpolate using Lagrange's formula */
	for( y_i = 0; y_i < 4; y_i++ ) {
		for( x_i = 0; x_i < 4; x_i++ ) {
			nodeIndex[x_i][y_i] = x_0 + x_i + ( y_0 + y_i ) * mesh->nVerts[0];
		}
	}

	for( x_i = 0; x_i < 4; x_i++ ) {
		px[x_i] = mesh->verts[nodeIndex[x_i][0]][0];
	}
	for( y_i = 0; y_i < 4; y_i++ ) {
		py[y_i] = mesh->verts[nodeIndex[0][y_i]][1];
	}
	
	for( y_i = 0; y_i < 4; y_i++ ) {
		for( x_i = 0; x_i < 4; x_i++ ) {
			for( dof_i = 0; dof_i < field->nDofs; dof_i++ ) {
				ptsX[x_i][dof_i] = field->vals[nodeIndex[x_i][y_i]][dof_i];
			}
		}

		InterpLagrange( coord[0], px, field->nDofs, ptsX, ptsY[y_i] );
	}

	InterpLagrange( coord[1], py, field->nDofs, ptsY, psi );
}

void Advector::ExtrapolateVelocities() {
	int node_i, dof_i;

	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		for( dof_i = 0; dof_i < velocity->nDofs; dof_i++ ) {
			velPlusOne->vals[node_i][dof_i]   = 2.0*velocity->vals[node_i][dof_i] - 1.0*velMinusOne->vals[node_i][dof_i];
			velPlusHalf->vals[node_i][dof_i]  = 1.5*velocity->vals[node_i][dof_i] - 0.5*velMinusOne->vals[node_i][dof_i];
			velMinusHalf->vals[node_i][dof_i] = 0.5*velocity->vals[node_i][dof_i] + 0.5*velMinusOne->vals[node_i][dof_i];
		}
	}
}
