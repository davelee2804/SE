#include <iostream>
#include <fstream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "Params.h"
#include "2LRLOps.h"

using namespace std;
using std::string;

BetaInvMatrix::BetaInvMatrix( string _name, Field* _rowField, Field* _colField, double _constant, double _alpha, double _f0, double _beta ) :
        Operator( _name, _rowField, _colField, _constant )
{
	alpha = _alpha;		/* BDF constant */
        f0    = _f0;
        beta  = _beta;
}

BetaInvMatrix::~BetaInvMatrix() {}

void BetaInvMatrix::AssembleElement( int el_i, double* M ) {
        int row, col, node_i;
        double *coord, weight, gCoord[3], detJac, f, *Ni, invFac, a00, a01, a10;
        Mesh* mesh = rowField->mesh;
        Element* el = mesh->el;
        int nNodes = el->nNodes;
        int nEntries = 2*nNodes;

        for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		coord = el->quadPts[pt_i]->coord;
		weight = el->quadPts[pt_i]->weight;
		detJac = mesh->DetJac( el_i, pt_i );
		Ni = el->ShapeFuncs( pt_i );
		mesh->LocalToGlobal( coord, el_i, gCoord );
		f = f0 + beta*gCoord[1];
		invFac = 1.0/(detJac*weight*detJac*weight*( alpha*alpha + constant*constant*f*f ));
		node_i = pt_i;
		row = node_i*2;
		col = node_i*2;
		a00 = +detJac*weight*alpha*Ni[node_i]*Ni[node_i];
		a01 = -detJac*weight*constant*f*Ni[node_i]*Ni[node_i];
		a10 = +detJac*weight*constant*f*Ni[node_i]*Ni[node_i];
		M[(row+0)*nEntries + (col+0)] += invFac*a00;
                M[(row+0)*nEntries + (col+1)] -= invFac*a01;
                M[(row+1)*nEntries + (col+0)] -= invFac*a10;
                M[(row+1)*nEntries + (col+1)] += invFac*a00;
        }
}

PhiMatrix::PhiMatrix( string _name, Field* _rowField, Field* _colField, double _constant, double _alpha, double _dt, double _f0, double _beta ) :
        Operator( _name, _rowField, _colField, _constant )
{
	alpha = _alpha;		/* BDF constant */
	dt    = _dt;
        f0    = _f0;
        beta  = _beta;
}

PhiMatrix::~PhiMatrix() {}

void PhiMatrix::AssembleElement( int el_i, double* M ) {
        double *coord, weight, gCoord[3], detJac, f, **GNix, invFac;
        Mesh* mesh = rowField->mesh;
        Element* el = mesh->el;
        int nNodes = el->nNodes;

        for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		coord 	= el->quadPts[pt_i]->coord;
		weight 	= el->quadPts[pt_i]->weight;
		GNix   	= mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		mesh->LocalToGlobal( coord, el_i, gCoord );
		f 	= f0 + beta*gCoord[1];
		invFac 	= 1.0/( alpha*alpha + dt*dt*f*f );
		for( int row_i = 0; row_i < el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < el->nNodes; col_j++ ) {
				M[row_i*nNodes + col_j] -= detJac*weight*constant*invFac*( 
					GNix[0][row_i]*( alpha*GNix[0][col_j] + dt*f*GNix[1][col_j] ) + 
					GNix[1][row_i]*( alpha*GNix[1][col_j] - dt*f*GNix[0][col_j] ) );
			}
		}
        }
}
/* note: doesn't include contribution from SVV (for now) - based on linear single layer viscosity */
GLayerVector::GLayerVector( string _name, Mesh* _mesh, double _constant, Field* _field, 
			    Field* _height1, Field* _height2, Field* _velocity1, Field* _velocity2, 
			    double _H1, double _H2, double _nu, double _gamma, double _gPrime, double _tau0, double _kws ) : 
RHSOp( _name, _mesh, _constant, _field ) 
{
	height1		= _height1;
	height2		= _height2;
	velocity1	= _velocity1;
	velocity2	= _velocity2;
	H1		= _H1;
	H2		= _H2;
	nu		= _nu;
	gamma		= _gamma;
	gPrime		= _gPrime;
	tau0		= _tau0;
	kws		= _kws;

	db     = new double*[1];  db[0]  = new double[2];
	dh1    = new double*[1];  dh1[0] = new double[2];
	dh2    = new double*[1];  dh2[0] = new double[2];
	du1    = new double*[2];  du1[0] = new double[2];  du1[1] = new double[2];
	du2    = new double*[2];  du2[0] = new double[2];  du2[1] = new double[2];

	db[0][0] = db[0][1] = 0.0;
}

GLayerVector::~GLayerVector() {
	delete[] db[0];   delete[] db;
	delete[] dh1[0];  delete[] dh1;
	delete[] dh2[0];  delete[] dh2;
	delete[] du1[0];  delete[] du1[1];  delete[] du1;
	delete[] du2[0];  delete[] du2[1];  delete[] du2;
}

void GLayerVector::AssembleElement( int el_i, double* G ) {
	double *coord, weight, detJac, *Ni, **GNix, h1, u1[2], h2, u2[2], a, b = 0.0, gCoord[2];

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		Ni     = mesh->el->ShapeFuncs( pt_i );
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );

		mesh->LocalToGlobal( coord, el_i, gCoord );

		if( field != NULL ) {
			field->InterpLocal( el_i, coord, &b );
			field->InterpDerivsGlobal( gCoord, db );
		}
		height1->InterpLocal( el_i, coord, &h1 );
		height2->InterpLocal( el_i, coord, &h2 );
		height1->InterpDerivsGlobal( gCoord, dh1 );
		height2->InterpDerivsGlobal( gCoord, dh2 );
		velocity1->InterpLocal( el_i, coord, u1 );
		velocity2->InterpLocal( el_i, coord, u2 );
		velocity1->InterpDerivsGlobal( gCoord, du1 );
		velocity2->InterpDerivsGlobal( gCoord, du2 );

		a = detJac*weight*constant/( H1 + H2 - b );

		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			/* non linear terms */
			//G[node_i*2+0] -= a*( GNix[0][node_i]*h1*u1[0]*u1[0] + GNix[1][node_i]*h1*u1[0]*u1[1] );
			//G[node_i*2+1] -= a*( GNix[1][node_i]*h1*u1[1]*u1[1] + GNix[0][node_i]*h1*u1[0]*u1[1] );
			//G[node_i*2+0] -= a*( GNix[0][node_i]*h2*u2[0]*u2[0] + GNix[1][node_i]*h2*u2[0]*u2[1] );
			//G[node_i*2+1] -= a*( GNix[1][node_i]*h2*u2[1]*u2[1] + GNix[0][node_i]*h2*u2[0]*u2[1] );
			G[node_i*2+0] += a*( 2.0*h1*u1[0]*du1[0][0] + u1[0]*u1[0]*dh1[0][0] + h1*u1[1]*du1[0][1] + u1[0]*u1[1]*dh1[0][1] + h1*u1[1]*du1[0][1] )*Ni[node_i];
			G[node_i*2+1] += a*( 2.0*h1*u1[1]*du1[1][1] + h1*u1[0]*du1[1][0] + u1[0]*u1[1]*dh1[0][0] + u1[1]*u1[1]*dh1[0][1] + h1*u1[0]*du1[1][0] )*Ni[node_i];
			G[node_i*2+0] += a*( 2.0*h2*u2[0]*du2[0][0] + u2[0]*u2[0]*dh2[0][0] + h2*u2[1]*du2[0][1] + u2[0]*u2[1]*dh2[0][1] + h2*u2[1]*du2[0][1] )*Ni[node_i];
			G[node_i*2+1] += a*( 2.0*h2*u2[1]*du2[1][1] + h2*u2[0]*du2[1][0] + u2[0]*u2[1]*dh2[0][0] + u2[1]*u2[1]*dh2[0][1] + h2*u2[0]*du2[1][0] )*Ni[node_i];


			/* viscous terms */
			G[node_i*2+0] += a*nu*h1*( 2.0*GNix[0][node_i]*du1[0][0] + GNix[1][node_i]*( du1[0][1] + du1[1][0] ) );
			G[node_i*2+1] += a*nu*h1*( 2.0*GNix[1][node_i]*du1[1][1] + GNix[0][node_i]*( du1[0][1] + du1[1][0] ) );
			G[node_i*2+0] += a*nu*h2*( 2.0*GNix[0][node_i]*du2[0][0] + GNix[1][node_i]*( du2[0][1] + du2[1][0] ) );
			G[node_i*2+1] += a*nu*h2*( 2.0*GNix[1][node_i]*du2[1][1] + GNix[0][node_i]*( du2[0][1] + du2[1][0] ) );

			/* bottom friction */
			G[node_i*2+0] += a*gamma*h2*u2[0]*Ni[node_i];
			G[node_i*2+1] += a*gamma*h2*u2[1]*Ni[node_i];

			/* bottom layer interface and topography gradient */
			G[node_i*2+0] += a*gPrime*h2*dh2[0][0]*Ni[node_i];
			G[node_i*2+1] += a*gPrime*h2*dh2[0][1]*Ni[node_i];
			//G[node_i*2+0] += a*gPrime*h2*db[0][0]*Ni[node_i];
			//G[node_i*2+1] += a*gPrime*h2*db[0][1]*Ni[node_i];

			/* top layer wind stress */
			//G[node_i*2+0] -= a*tau0*cos( kws*gCoord[1] )*Ni[node_i];
			G[node_i*2+0] += a*tau0*sin( kws*gCoord[1] )*Ni[node_i];
		}
	}
}

GLayerVector2::GLayerVector2( string _name, Mesh* _mesh, double _constant, Field* _field, Field* _velocity1, Field* _velocity2, 
			    double _H1, double _H2, double _nu, double _gPrime, double _tau0, double _kws ) : 
RHSOp( _name, _mesh, _constant, _field ) 
{
	velocity1	= _velocity1;
	velocity2	= _velocity2;
	H1		= _H1;
	H2		= _H2;
	nu		= _nu;
	gPrime		= _gPrime;
	tau0		= _tau0;
	kws		= _kws;

	dei    = new double*[1];  dei[0] = new double[2];
	du1    = new double*[2];  du1[0] = new double[2];  du1[1] = new double[2];
	du2    = new double*[2];  du2[0] = new double[2];  du2[1] = new double[2];
}

GLayerVector2::~GLayerVector2() {
	delete[] dei[0];  delete[] dei;
	delete[] du1[0];  delete[] du1[1];  delete[] du1;
	delete[] du2[0];  delete[] du2[1];  delete[] du2;
}

void GLayerVector2::AssembleElement( int el_i, double* G ) {
	double *coord, weight, detJac, *Ni, **GNix, h1, u1[2], h2, u2[2], a, e, gCoord[2];

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		Ni     = mesh->el->ShapeFuncs( pt_i );
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );

		mesh->LocalToGlobal( coord, el_i, gCoord );

		field->InterpLocal( el_i, coord, &e );
		field->InterpDerivsGlobal( gCoord, dei );
		velocity1->InterpLocal( el_i, coord, u1 );
		velocity2->InterpLocal( el_i, coord, u2 );
		velocity1->InterpDerivsGlobal( gCoord, du1 );
		velocity2->InterpDerivsGlobal( gCoord, du2 );

		h1 = H1 - e;
		h2 = H2 + e;

		a = detJac*weight*constant/( H1 + H2 );

		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			/* non linear terms */
			//G[node_i*2+0] += a*( 2.0*h1*u1[0]*du1[0][0] + u1[0]*u1[0]*dh1[0][0] + h1*u1[1]*du1[0][1] + u1[0]*u1[1]*dh1[0][1] + h1*u1[1]*du1[0][1] )*Ni[node_i];
			//G[node_i*2+1] += a*( 2.0*h1*u1[1]*du1[1][1] + h1*u1[0]*du1[1][0] + u1[0]*u1[1]*dh1[0][0] + u1[1]*u1[1]*dh1[0][1] + h1*u1[0]*du1[1][0] )*Ni[node_i];
			//G[node_i*2+0] += a*( 2.0*h2*u2[0]*du2[0][0] + u2[0]*u2[0]*dh2[0][0] + h2*u2[1]*du2[0][1] + u2[0]*u2[1]*dh2[0][1] + h2*u2[1]*du2[0][1] )*Ni[node_i];
			//G[node_i*2+1] += a*( 2.0*h2*u2[1]*du2[1][1] + h2*u2[0]*du2[1][0] + u2[0]*u2[1]*dh2[0][0] + u2[1]*u2[1]*dh2[0][1] + h2*u2[0]*du2[1][0] )*Ni[node_i];
			G[node_i*2+0] += a*( 2.0*h1*u1[0]*du1[0][0] - u1[0]*u1[0]*dei[0][0] + h1*u1[1]*du1[0][1] - u1[0]*u1[1]*dei[0][1] + h1*u1[1]*du1[0][1] )*Ni[node_i];
			G[node_i*2+1] += a*( 2.0*h1*u1[1]*du1[1][1] + h1*u1[0]*du1[1][0] - u1[0]*u1[1]*dei[0][0] - u1[1]*u1[1]*dei[0][1] + h1*u1[0]*du1[1][0] )*Ni[node_i];
			G[node_i*2+0] += a*( 2.0*h2*u2[0]*du2[0][0] + u2[0]*u2[0]*dei[0][0] + h2*u2[1]*du2[0][1] + u2[0]*u2[1]*dei[0][1] + h2*u2[1]*du2[0][1] )*Ni[node_i];
			G[node_i*2+1] += a*( 2.0*h2*u2[1]*du2[1][1] + h2*u2[0]*du2[1][0] + u2[0]*u2[1]*dei[0][0] + u2[1]*u2[1]*dei[0][1] + h2*u2[0]*du2[1][0] )*Ni[node_i];


			/* viscous terms */
			G[node_i*2+0] += a*nu*h1*( 2.0*GNix[0][node_i]*du1[0][0] + GNix[1][node_i]*( du1[0][1] + du1[1][0] ) );
			G[node_i*2+1] += a*nu*h1*( 2.0*GNix[1][node_i]*du1[1][1] + GNix[0][node_i]*( du1[0][1] + du1[1][0] ) );
			G[node_i*2+0] += a*nu*h2*( 2.0*GNix[0][node_i]*du2[0][0] + GNix[1][node_i]*( du2[0][1] + du2[1][0] ) );
			G[node_i*2+1] += a*nu*h2*( 2.0*GNix[1][node_i]*du2[1][1] + GNix[0][node_i]*( du2[0][1] + du2[1][0] ) );

			/* interface gradient */
			G[node_i*2+0] += a*gPrime*h2*dei[0][0]*Ni[node_i];
			G[node_i*2+1] += a*gPrime*h2*dei[0][1]*Ni[node_i];

			/* top layer wind stress */
			G[node_i*2+0] += a*tau0*sin( kws*gCoord[1] )*Ni[node_i];
		}
	}
}

DivTotHeightVelRHS::DivTotHeightVelRHS( string _name, Mesh* _mesh, double _constant, Field* _field, Field* _velocity, double _H, Field* _topo ) :
RHSOp( _name, _mesh, _constant, _field ) 
{
	velocity 	= _velocity;
	H		= _H;
	topo		= _topo;
}

DivTotHeightVelRHS::~DivTotHeightVelRHS() {}

void DivTotHeightVelRHS::AssembleElement( int el_i, double* rhs ) {
	double	*coord, weight, detJac, **GNix, h, t, v[2], a;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpLocal( el_i, coord, &h );
		velocity->InterpLocal( el_i, coord, v );
		topo->InterpLocal( el_i, coord, &t );
		a = detJac*weight*constant*( H + h - t );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] -= a*( GNix[0][node_i]*v[0] + GNix[1][node_i]*v[1] );
		}
	}
}

BetaRHS::BetaRHS( string _name, Mesh* _mesh, double _constant, Field* _field, double _f0, double _beta ) : RHSOp( _name, _mesh, _constant, _field ) {
	f0 	= _f0;
	beta	= _beta;
}

BetaRHS::~BetaRHS() {}

void BetaRHS::AssembleElement( int el_i, double* rhs ) {
	double detJac, weight, *coord, *Ni, u[2], gCoord[2], a;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		weight = mesh->el->quadPts[pt_i]->weight;
		coord  = mesh->el->quadPts[pt_i]->coord;
		detJac = mesh->DetJac( el_i, pt_i );
		Ni     = mesh->el->ShapeFuncs( pt_i );
		field->InterpLocal( el_i, coord, u );
		mesh->LocalToGlobal( coord, el_i, gCoord );
		a      = detJac*weight*constant*( f0 + beta*gCoord[1] );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[2*node_i+0] -= a*u[1]*Ni[node_i];
			rhs[2*node_i+1] += a*u[0]*Ni[node_i];
		}
	}
}

BaroclinicPressureBC::BaroclinicPressureBC( string _name, Mesh* _mesh, double _constant, Field* _field, Params* _p, 
					    Field* _velTop, Field* _velBot, Field* _velBar, Field* _height1, Field* _height2 ) : 
					    RHSOp( _name, _mesh, _constant, _field ) 
{
	p = _p;
	velTop = _velTop;
	velBot = _velBot;
	velBar = _velBar;
	height1 = _height1;
	height2 = _height2;

	dh = new double*[1];
	dh[0] = new double[2];
	du1 = new double*[2];
	du1[0] = new double[2];
	du1[1] = new double[2];
	du2 = new double*[2];
	du2[0] = new double[2];
	du2[1] = new double[2];
}

BaroclinicPressureBC::~BaroclinicPressureBC() {
	delete[] dh[0];
	delete[] dh;
	delete[] du1[0];
	delete[] du1[1];
	delete[] du1;
	delete[] du2[0];
	delete[] du2[1];
	delete[] du2;
}

void BaroclinicPressureBC::AssembleBottom( int el_i, double* rhs ) {
	double	detJac	= fabs( mesh->max[0] - mesh->min[0] )/(2.0*mesh->nEls[0]);
	//double	x[2];
	double	norm	= -1.0;
	double	coord[2];
	double	gCoord[2];
	//double	tau;
	//double	h;
	double	weight;
	double	Ni	= 1.0;	/* shape functions are orthoginal to quadrature points for legendre basis */
	int	start	= 0;
	int	skip	= 1;
	double	v[2];
	double	f;
	double	u1[2], u2[2];
	double	h1, h2;

	//x[1] = mesh->min[1];
	coord[1] = -1.0;

	for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
		coord[0] = mesh->el->abcissa[pt_i];
		weight = mesh->el->weights[pt_i];
		mesh->LocalToGlobal( coord, el_i, gCoord );
		height1->InterpLocal( el_i, coord, &h1 );
		height2->InterpLocal( el_i, coord, &h2 );
		height2->InterpDerivsGlobal( gCoord, dh ); //TODO: add topography gradients
		//tau = 0.0;
		velTop->InterpLocal( el_i, coord, u1 );
		velBot->InterpLocal( el_i, coord, u2 );
		velTop->InterpDerivsGlobal( gCoord, du1 );
		velBot->InterpDerivsGlobal( gCoord, du2 );
		velBar->InterpLocal( el_i, coord, v );
		f = p->H*( p->f0 + p->beta*gCoord[1] );
		rhs[start + skip*pt_i] += constant*norm*detJac*weight*Ni*( /*tau +*/ h2*dh[0][1] + f*v[0] + h1*u1[0]*du1[1][0] + h2*u2[0]*du2[1][0] );
	}
}

void BaroclinicPressureBC::AssembleTop( int el_i, double* rhs ) {
	double	detJac	= fabs( mesh->max[0] - mesh->min[0] )/(2.0*mesh->nEls[0]);
	//double	x[2];
	double	norm	= +1.0;
	double	coord[2];
	double	gCoord[2];
	//double	tau;
	//double	h;
	double	weight;
	double	Ni	= 1.0;	/* shape functions are orthoginal to quadrature points for legendre basis */
	int	start	= (mesh->el->N+1)*(mesh->el->N);
	int	skip	= 1;
	double	v[2];
	double	f;
	double	u1[2], u2[2];
	double	h1, h2;

	//x[1] = mesh->max[1];
	coord[1] = +1.0;

	for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
		coord[0] = mesh->el->abcissa[pt_i];
		weight = mesh->el->weights[pt_i];
		mesh->LocalToGlobal( coord, el_i, gCoord );
		height1->InterpLocal( el_i, coord, &h1 );
		height2->InterpLocal( el_i, coord, &h2 );
		height2->InterpDerivsGlobal( gCoord, dh ); //TODO: add topography gradients
		//tau = 0.0;
		velTop->InterpLocal( el_i, coord, u1 );
		velBot->InterpLocal( el_i, coord, u2 );
		velTop->InterpDerivsGlobal( gCoord, du1 );
		velBot->InterpDerivsGlobal( gCoord, du2 );
		velBar->InterpLocal( el_i, coord, v );
		f = p->H*( p->f0 + p->beta*gCoord[1] );
		rhs[start + skip*pt_i] += constant*norm*detJac*weight*Ni*( /*tau +*/ h2*dh[0][1] + f*v[0] + h1*u1[0]*du1[1][0] + h2*u2[0]*du2[1][0] );
	}
}

void BaroclinicPressureBC::AssembleLeft( int el_i, double* rhs ) {
	double	detJac	= fabs( mesh->max[1] - mesh->min[1] )/(2.0*mesh->nEls[1]);
	//double	x[2];
	double	norm	= -1.0;
	double	coord[2];
	double	gCoord[2];
	//double	tau;
	//double	h;
	double	weight;
	double	Ni	= 1.0;	/* shape functions are orthoginal to quadrature points for legendre basis */
	int	start	= 0;
	int	skip	= mesh->el->N+1;
	double	v[2];
	double	f;
	double	u1[2], u2[2];
	double	h1, h2;

	//x[0] = mesh->min[0];
	coord[0] = -1.0;

	for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
		coord[1] = mesh->el->abcissa[pt_i];
		weight = mesh->el->weights[pt_i];
		mesh->LocalToGlobal( coord, el_i, gCoord );
		height1->InterpLocal( el_i, coord, &h1 );
		height2->InterpLocal( el_i, coord, &h2 );
		height2->InterpDerivsGlobal( gCoord, dh ); //TODO: add topography gradients
		//tau = p->tau*cos( p->k*gCoord[1] );
		velTop->InterpLocal( el_i, coord, u1 );
		velBot->InterpLocal( el_i, coord, u2 );
		velTop->InterpDerivsGlobal( gCoord, du1 );
		velBot->InterpDerivsGlobal( gCoord, du2 );
		velBar->InterpLocal( el_i, coord, v );
		f = p->H*( p->f0 + p->beta*gCoord[1] );
		rhs[start + skip*pt_i] += constant*norm*detJac*weight*Ni*( /*tau +*/ h2*dh[0][0] - f*v[1] + h1*u1[1]*du1[0][1] + h2*u2[1]*du2[0][1] );
	}
}

void BaroclinicPressureBC::AssembleRight( int el_i, double* rhs ) {
	double	detJac	= fabs( mesh->max[1] - mesh->min[1] )/(2.0*mesh->nEls[1]);
	//double	x[2];
	double	norm	= +1.0;
	double	coord[2];
	double	gCoord[2];
	//double	tau;
	//double	h;
	double	weight;
	double	Ni	= 1.0;	/* shape functions are orthoginal to quadrature points for legendre basis */
	int	start	= mesh->el->N;
	int	skip	= mesh->el->N+1;
	double	v[2];
	double	f;
	double	u1[2], u2[2];
	double	h1, h2;

	//x[0] = mesh->max[0];
	coord[0] = +1.0;

	for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
		coord[1] = mesh->el->abcissa[pt_i];
		weight = mesh->el->weights[pt_i];
		mesh->LocalToGlobal( coord, el_i, gCoord );
		height1->InterpLocal( el_i, coord, &h1 );
		height2->InterpLocal( el_i, coord, &h2 );
		height2->InterpDerivsGlobal( gCoord, dh ); //TODO: add topography gradients
		//tau = p->tau*cos( p->k*gCoord[1] );
		velTop->InterpLocal( el_i, coord, u1 );
		velBot->InterpLocal( el_i, coord, u2 );
		velTop->InterpDerivsGlobal( gCoord, du1 );
		velBot->InterpDerivsGlobal( gCoord, du2 );
		velBar->InterpLocal( el_i, coord, v );
		f = p->H*( p->f0 + p->beta*gCoord[1] );
		rhs[start + skip*pt_i] += constant*norm*detJac*weight*Ni*( /*tau +*/ h2*dh[0][0] - f*v[1] + h1*u1[1]*du1[0][1] + h2*u2[1]*du2[0][1] );
	}
}

void BaroclinicPressureBC::AssembleElement( int el_i, double* rhs ) {
	if( el_i/mesh->nEls[0] == 0 )                 { AssembleBottom( el_i, rhs ); }
	if( el_i/mesh->nEls[0] == mesh->nEls[1] - 1 ) { AssembleTop( el_i, rhs ); }
	if( el_i%mesh->nEls[0] == 0 )                 { AssembleLeft( el_i, rhs ); }
	if( el_i%mesh->nEls[0] == mesh->nEls[0] - 1 ) { AssembleRight( el_i, rhs ); }
}

BarotropicPressureBC::BarotropicPressureBC( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {}

BarotropicPressureBC::~BarotropicPressureBC() {}

void BarotropicPressureBC::AssembleElement( int el_i, double* rhs ) {
	if( el_i/mesh->nEls[0] == 0 )                 { AssembleBottom( el_i, rhs ); }
	if( el_i/mesh->nEls[0] == mesh->nEls[1] - 1 ) { AssembleTop( el_i, rhs ); }
	if( el_i%mesh->nEls[0] == 0 )                 { AssembleLeft( el_i, rhs ); }
	if( el_i%mesh->nEls[0] == mesh->nEls[0] - 1 ) { AssembleRight( el_i, rhs ); }
}

void BarotropicPressureBC::AssembleBottom( int el_i, double* rhs ) {
	double	detJac	= fabs( mesh->max[0] - mesh->min[0] )/(2.0*mesh->nEls[0]);
	//double	x[2];
	double	norm	= -1.0;
	double	coord[2];
	double	weight;
	double	Ni	= 1.0;	/* shape functions are orthoginal to quadrature points for legendre basis */
	int	start	= 0;
	int	skip	= 1;
	double	v[2];

	//x[1] = mesh->min[1];
	coord[1] = -1.0;

	for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
		coord[0] = mesh->el->abcissa[pt_i];
		weight = mesh->el->weights[pt_i];
		field->InterpLocal( el_i, coord, v );
		rhs[start + skip*pt_i] += constant*norm*detJac*weight*v[1]*Ni;
	}
}

void BarotropicPressureBC::AssembleTop( int el_i, double* rhs ) {
	double	detJac	= fabs( mesh->max[0] - mesh->min[0] )/(2.0*mesh->nEls[0]);
	//double	x[2];
	double	norm	= +1.0;
	double	coord[2];
	double	weight;
	double	Ni	= 1.0;	/* shape functions are orthoginal to quadrature points for legendre basis */
	int	start	= (mesh->el->N+1)*(mesh->el->N);
	int	skip	= 1;
	double	v[2];

	//x[1] = mesh->max[1];
	coord[1] = +1.0;

	for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
		coord[0] = mesh->el->abcissa[pt_i];
		weight = mesh->el->weights[pt_i];
		field->InterpLocal( el_i, coord, v );
		rhs[start + skip*pt_i] += constant*norm*detJac*weight*v[1]*Ni;
	}
}

void BarotropicPressureBC::AssembleLeft( int el_i, double* rhs ) {
	double	detJac	= fabs( mesh->max[1] - mesh->min[1] )/(2.0*mesh->nEls[1]);
	//double	x[2];
	double	norm	= -1.0;
	double	coord[2];
	double	weight;
	double	Ni	= 1.0;	/* shape functions are orthoginal to quadrature points for legendre basis */
	int	start	= 0;
	int	skip	= mesh->el->N+1;
	double	v[2];

	//x[0] = mesh->min[0];
	coord[0] = -1.0;

	for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
		coord[1] = mesh->el->abcissa[pt_i];
		weight = mesh->el->weights[pt_i];
		field->InterpLocal( el_i, coord, v );
		rhs[start + skip*pt_i] += constant*norm*detJac*weight*v[0]*Ni;
	}
}

void BarotropicPressureBC::AssembleRight( int el_i, double* rhs ) {
	double	detJac	= fabs( mesh->max[1] - mesh->min[1] )/(2.0*mesh->nEls[1]);
	//double	x[2];
	double	norm	= +1.0;
	double	coord[2];
	double	weight;
	double	Ni	= 1.0;	/* shape functions are orthoginal to quadrature points for legendre basis */
	int	start	= mesh->el->N;
	int	skip	= mesh->el->N+1;
	double	v[2];

	//x[0] = mesh->max[0];
	coord[0] = +1.0;

	for( int pt_i = 0; pt_i < mesh->el->N+1; pt_i++ ) {
		coord[1] = mesh->el->abcissa[pt_i];
		weight = mesh->el->weights[pt_i];
		field->InterpLocal( el_i, coord, v );
		rhs[start + skip*pt_i] += constant*norm*detJac*weight*v[0]*Ni;
	}
}

WindStress2RHS::WindStress2RHS( string _name, Mesh* _mesh, double _constant, Field* _field, double _k, double _H ) : RHSOp( _name, _mesh, _constant, _field ) {
	k = _k;
	H = _H;
}

WindStress2RHS::~WindStress2RHS() {}

void WindStress2RHS::AssembleElement( int el_i, double* rhs ) {
	double	*Ni, e = 0.0, *coord, detJac, weight, gCoord[2];

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		detJac = mesh->DetJac( el_i, pt_i );
		Ni     = mesh->el->ShapeFuncs( pt_i );
		mesh->LocalToGlobal( coord, el_i, gCoord );
		if( field ) {
			field->InterpLocal( el_i, coord, &e );
		}
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[2*node_i+0] += detJac*weight*( constant*sin( k*gCoord[1] )/( H - e ) )*Ni[node_i];
		}
	}
}

PhiRHS::PhiRHS( string _name, Mesh* _mesh, double _constant, Field* _field, 
		double _sign, double _dt, double _H, double _tau0, double _k, double _nu, double _p, double _f0, double _beta,
		Field* _pres, Field* _eta, Field* _fieldSL ) : 
RHSOp( _name, _mesh, _constant, _field ) {
	sign	= _sign;
	dt	= _dt;
	H 	= _H;
	tau0 	= _tau0;
	k 	= _k;
	nu 	= _nu;
	p	= _p;
	f0	= _f0;
	beta	= _beta;
	pres 	= _pres;
	eta 	= _eta;
	fieldSL = _fieldSL;

	d2u    = new double*[2];
	d2u[0] = new double[3];
	d2u[1] = new double[3];
	dp     = new double*[1];
	dp[0]  = new double[2];
}

PhiRHS::~PhiRHS() {
	delete[] d2u[0];
	delete[] d2u[1];
	delete[] d2u;
	delete[] dp[0];
	delete[] dp;
}

void PhiRHS::AssembleElement( int el_i, double* rhs ) {
	double *coord, weight, detJac, **GNix;
	double usl[2], e, gCoord[2], f, detInv, a[2];

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		mesh->LocalToGlobal( coord, el_i, gCoord );
		f      = f0 + beta*gCoord[1];
		detInv = 1.0/(constant*constant + dt*dt*f*f);
		fieldSL->InterpLocal( el_i, coord, usl );
		pres->InterpDerivsGlobal( gCoord, dp );
		eta->InterpLocal( el_i, coord, &e );
		field->InterpSecondDerivsAtCoord( el_i, coord, d2u ); //TODO: test this func!!
		a[0] = usl[0] - dt*p*dp[0][0] + dt*nu*d2u[0][0] + dt*tau0*sin(k*gCoord[1])/(H-e);
		a[1] = usl[1] - dt*p*dp[0][1] + dt*nu*d2u[1][1];
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] -= 0.5*detJac*weight*sign*dt*H*detInv*GNix[0][node_i]*(constant*a[0] + dt*f*a[1]);
			rhs[node_i] -= 0.5*detJac*weight*sign*dt*H*detInv*GNix[1][node_i]*(constant*a[1] - dt*f*a[0]);
		}
	}
}
