#include <cstdlib>
#include <string>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "QuadPoint.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Operator.h"
#include "RHSOp.h"
#include "SWOps.h"

using namespace std;
using std::string;

DivPhiVelMatrix::DivPhiVelMatrix( string _name, Field* _rowField, Field* _colField, double _constant, Field* _field ) : Operator( _name, _rowField, _colField, _constant ) {
	field = _field;
}

DivPhiVelMatrix::~DivPhiVelMatrix() {}

void DivPhiVelMatrix::AssembleElement( int el_i, double* M ) {
	int pt_i, row_i, col_j;
	double *coord, weight, detJac, *Ni, **GNix, v[2];
	Mesh* mesh = rowField->mesh;
	Element* el = mesh->el;

	for( pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		coord  = el->quadPts[pt_i]->coord;
		weight = el->quadPts[pt_i]->weight;
		Ni     = el->ShapeFuncs( pt_i );
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpLocal( el_i, coord, v );
		for( row_i = 0; row_i < el->nNodes; row_i++ ) {
			for( col_j = 0; col_j < el->nNodes; col_j++ ) {
				M[row_i*el->nNodes+col_j] += detJac*weight*( GNix[0][row_i]*v[0] + GNix[1][row_i]*v[1] )*Ni[col_j];
			}
		}
	}
}

WindStressRHS::WindStressRHS( string _name, Mesh* _mesh, double _constant, Field* _field, double _k, double _H ) : RHSOp( _name, _mesh, _constant, _field ) {
	k = _k;
	H = _H;
}

WindStressRHS::~WindStressRHS() {}

void WindStressRHS::AssembleElement( int el_i, double* rhs ) {
	double	*Ni, h, *coord, detJac, weight, gCoord[2];

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		detJac = mesh->DetJac( el_i, pt_i );
		Ni     = mesh->el->ShapeFuncs( pt_i );
		mesh->LocalToGlobal( coord, el_i, gCoord );
		field->InterpLocal( el_i, coord, &h );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[2*node_i+0] -= detJac*weight*( constant*cos( k*gCoord[1] )/( H - h ) )*Ni[node_i];
		}
	}
}

DivHeightMinusTopo::DivHeightMinusTopo( string _name, Field* _rowField, Field* _colField, double _constant, double _H, Field* _topo ) : Operator( _name, _rowField, _colField, _constant ) {
	H = _H;
	topo = _topo;

	db = new double*[1];
	db[0] = new double[2];
}

DivHeightMinusTopo::~DivHeightMinusTopo() {
	delete[] db[0];
	delete[] db;
}

void DivHeightMinusTopo::AssembleElement( int el_i, double* M ) {
	double 	detJac, weight, *coord, *Ni, *Nj, **GNjx, b, a, gCoord[2];
	int	nCols 	= 2*colField->mesh->el->nNodes;

	for( int pt_i = 0; pt_i < rowField->mesh->el->nPoints; pt_i++ ) {
		coord  = rowField->mesh->el->quadPts[pt_i]->coord;
		weight = rowField->mesh->el->quadPts[pt_i]->weight;
		Ni     = rowField->mesh->el->ShapeFuncs( pt_i );
		Nj     = colField->mesh->el->ShapeFuncsAtCoord( coord );
		GNjx   = colField->mesh->ShapeFuncDerivsAtCoord( el_i, coord, &detJac );
		rowField->mesh->LocalToGlobal( coord, el_i, gCoord );
		topo->InterpLocal( el_i, coord, &b );
		topo->InterpDerivsGlobal( gCoord, db );
		a      = detJac*weight*constant;
		for( int row_i = 0; row_i < rowField->mesh->el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < colField->mesh->el->nNodes; col_j++ ) {
					M[row_i*nCols+(2*col_j+0)] += a*Ni[row_i]*( ( H - b )*GNjx[0][col_j] - db[0][0]*Nj[col_j] );
					M[row_i*nCols+(2*col_j+1)] += a*Ni[row_i]*( ( H - b )*GNjx[1][col_j] - db[0][1]*Nj[col_j] );
			}
		}
	}
}

DivVelBarMatrix::DivVelBarMatrix( string _name, Field* _rowField, Field* _colField, double _constant, Field* _velBar ) : Operator( _name, _rowField, _colField, _constant ) {
	velBar = _velBar;
}

DivVelBarMatrix::~DivVelBarMatrix() {}

void DivVelBarMatrix::AssembleElement( int el_i, double* M ) {
	double *coord, detJac, weight, *Ni, **GNix, v[2];

	for( int pt_i = 0; pt_i < rowField->mesh->el->nPoints; pt_i++ ) {
		coord  = rowField->mesh->el->quadPts[pt_i]->coord;
		weight = rowField->mesh->el->quadPts[pt_i]->weight;
		Ni     = rowField->mesh->el->ShapeFuncs( pt_i );
		GNix   = rowField->mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		velBar->InterpLocal( el_i, coord, v );
		for( int row_i = 0; row_i < rowField->mesh->el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < colField->mesh->el->nNodes; col_j++ ) {
				M[row_i*colField->mesh->el->nNodes + col_j] -= constant*detJac*weight*( GNix[row_i][0]*v[0] + GNix[row_i][1]*v[1] )*Ni[col_j];
			}
		}
	}
}

DivVelBarRHS::DivVelBarRHS( string _name, Mesh* _mesh, double _constant, Field* _field, Field* _topo, double _H_i ) : RHSOp( _name, _mesh, _constant, _field ) {
	topo   = _topo;
	H_i    = _H_i;
}

DivVelBarRHS::~DivVelBarRHS() {}

void DivVelBarRHS::AssembleElement( int el_i, double* rhs ) {
	double *coord, detJac, weight, **GNix, v[2], b;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpLocal( el_i, coord, v );
		topo->InterpLocal( el_i, coord, &b );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] += constant*detJac*weight*( H_i - b )*( GNix[node_i][0]*v[0] + GNix[node_i][1]*v[1] );
		}
	}
}
