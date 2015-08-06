#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>

#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Legendre.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Operator.h"

using namespace std;
using std::string;

Operator::Operator( std::string _name, Field* _rowField, Field* _colField, double _constant ) {
	name     = _name;
	rowField = _rowField;
	colField = _colField;
	constant = _constant;
}

Operator::~Operator() {
}

Biharmonic::Biharmonic( string _name, Field* _rowField, Field* _colField, double _constant ) : Operator( _name, _rowField, _colField, _constant ) {}

Biharmonic::~Biharmonic() {}

void Biharmonic::AssembleElement( int el_i, double* M ) {
	Mesh*		mesh		= rowField->mesh;
	Legendre*	el		= (Legendre*)mesh->el;
	double 		**GNixx, detJac, weight;

	/* internal element assembly */
	for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		weight = el->quadPts[pt_i]->weight;
		GNixx  = mesh->ShapeFuncSecondDerivs( el_i, pt_i, &detJac );
		for( int row_i = 0; row_i < el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < el->nNodes; col_j++ ) {
				M[row_i*el->nNodes+col_j] += detJac*weight*constant*( GNixx[0][row_i]*GNixx[0][col_j] + 
										      GNixx[0][row_i]*GNixx[1][col_j] +
										      GNixx[1][row_i]*GNixx[0][col_j] +
										      GNixx[1][row_i]*GNixx[1][col_j] );
			}
		}
	}
}

Laplacian::Laplacian( string _name, Field* _rowField, Field* _colField, double _constant ) : Operator( _name, _rowField, _colField, _constant ) {
	if( rowField != colField ) {
		cerr << "ERROR: attemting to instantiate a Laplacian operator for which row field != column field.\n";
		exit(1);
	}
}

Laplacian::~Laplacian() {}

void Laplacian::AssembleElement( int el_i, double* M ) {
	Mesh*		mesh	= rowField->mesh;
	Element*	el	= mesh->el;
	double		**GNix, detJac, weight, w;
	int		row, col, nDofs = rowField->nDofs, nEntries = nDofs*el->nNodes;

	for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		weight = el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		detJac = 1.0;
		w      = detJac*weight*constant;
		for( int row_i = 0; row_i < el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < el->nNodes; col_j++ ) {
				row = nDofs*row_i;
				col = nDofs*col_j;
				for( int dof = 0; dof < nDofs; dof++ ) {
					M[(row+dof)*nEntries + (col+dof)] += w*( GNix[0][row_i]*GNix[0][col_j] + GNix[1][row_i]*GNix[1][col_j] );
				}
			}
		}
	}
}

NonLinearLaplacian::NonLinearLaplacian( string _name, Field* _rowField, Field* _colField, double _constant, Field *_field ) : Operator( _name, _rowField, _colField, _constant ) {
	field = _field;
}

NonLinearLaplacian::~NonLinearLaplacian() {}

void NonLinearLaplacian::AssembleElement( int el_i, double* M ) {
        Mesh*           mesh    = rowField->mesh;
        Element*        el      = mesh->el;
        double          **GNix, detJac, weight, *coord, w, a;
        int             row, col, nDofs = rowField->nDofs, nEntries = nDofs*el->nNodes;

        for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		coord  = el->quadPts[pt_i]->coord;
                weight = el->quadPts[pt_i]->weight;
                GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpLocal( el_i, coord, &a );
                w      = detJac*weight*constant*a;
                for( int row_i = 0; row_i < el->nNodes; row_i++ ) {
                        for( int col_j = 0; col_j < el->nNodes; col_j++ ) {
                                row = nDofs*row_i;
                                col = nDofs*col_j;
                                for( int dof = 0; dof < nDofs; dof++ ) {
                                        M[(row+dof)*nEntries + (col+dof)] += w*( GNix[0][row_i]*GNix[0][col_j] +
                                                                                 GNix[1][row_i]*GNix[1][col_j] );
                                }
                        }
                }
        }
}

#define NDOFS 2

StressTensor::StressTensor( string _name, Field* _rowField, Field* _colField, double _constant ) : Operator( _name, _rowField, _colField, _constant ) {
	if( rowField != colField ) {
		cerr << "ERROR: attemting to instantiate a stress tensor operator for which row field != column field.\n";
		exit(1);
	}
}

StressTensor::~StressTensor() {}

/* ASSUMPTION! row field = column field, see Hughes, pages 79-89 for details */
void StressTensor::AssembleElement( int el_i, double* M ) {
	Mesh* 		mesh 	= rowField->mesh;
	Element* 	el 	= mesh->el;
	double 		**GNix, detJac, weight, w;
	int		row, col, nEntries = NDOFS*el->nNodes;

	for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		weight = el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		w      = detJac*weight*constant;
		for( int row_i = 0; row_i < el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < el->nNodes; col_j++ ) {
				row = NDOFS*row_i;
				col = NDOFS*col_j;
				M[(row+0)*nEntries + (col+0)] += w*( 2*GNix[0][row_i]*GNix[0][col_j] + GNix[1][row_i]*GNix[1][col_j] );
				M[(row+0)*nEntries + (col+1)] += w*GNix[1][row_i]*GNix[0][col_j];
				M[(row+1)*nEntries + (col+0)] += w*GNix[0][row_i]*GNix[1][col_j];
				M[(row+1)*nEntries + (col+1)] += w*( 2*GNix[1][row_i]*GNix[1][col_j] + GNix[0][row_i]*GNix[0][col_j] );
			}
		}
	}
}

MassMatrix::MassMatrix( string _name, Field* _rowField, Field* _colField, double _constant ) : Operator( _name, _rowField, _colField, _constant ) {}

MassMatrix::~MassMatrix() {}

/* ASSUMPTION! row vector = column vector */
void MassMatrix::AssembleElement( int el_i, double* M ) {
	Mesh*		mesh 	= rowField->mesh;
	Element*	el	= mesh->el;
	double		*Ni, detJac, weight;
	int		row, col, nDofs = rowField->nDofs, nEntries = nDofs*el->nNodes;

	for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		weight = el->quadPts[pt_i]->weight;
		Ni     = el->ShapeFuncs( pt_i );
		detJac = mesh->DetJac( el_i, pt_i );
		for( int row_i = 0; row_i < el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < el->nNodes; col_j++ ) {
				row = row_i*nDofs;
				col = col_j*nDofs;
				for( int dof = 0; dof < nDofs; dof++ ) {
					M[(row+dof)*nEntries + (col+dof)] += detJac*weight*constant*Ni[row_i]*Ni[col_j];
				}
			}
		}
	}
}

Gradient::Gradient( string _name, Field* _rowField, Field* _colField, double _constant ) : Operator( _name, _rowField, _colField, _constant ) {}

Gradient::~Gradient() {}

void Gradient::AssembleElement( int el_i, double* M ) {
	double *coord, detJac, weight, **GNix, *Ni;

	for( int pt_i = 0; pt_i < rowField->mesh->el->nPoints; pt_i++ ) {
		coord  = rowField->mesh->el->quadPts[pt_i]->coord;
		weight = rowField->mesh->el->quadPts[pt_i]->weight;
		GNix   = rowField->mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		Ni     = colField->mesh->el->ShapeFuncsAtCoord( coord );
		for( int row_i = 0; row_i < rowField->mesh->el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < colField->mesh->el->nNodes; col_j++ ) {
				M[(2*row_i+0)*colField->mesh->el->nNodes + col_j] -= detJac*weight*constant*GNix[0][row_i]*Ni[col_j];
				M[(2*row_i+1)*colField->mesh->el->nNodes + col_j] -= detJac*weight*constant*GNix[1][row_i]*Ni[col_j];
			}
		}
	}
}

Divergence::Divergence( string _name, Field* _rowField, Field* _colField, double _constant ) : Operator( _name, _rowField, _colField, _constant ) {}

Divergence::~Divergence() {}

void Divergence::AssembleElement( int el_i, double* M ) {
/* this works, but poorer convergence... */
#if 0
	double *coord, detJac, weight, **GNix, *Ni;

	for( int pt_i = 0; pt_i < colField->mesh->el->nPoints; pt_i++ ) {
		coord  = colField->mesh->el->quadPts[pt_i]->coord;
		weight = colField->mesh->el->quadPts[pt_i]->weight;
		GNix   = colField->mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		Ni     = rowField->mesh->el->ShapeFuncsAtCoord( coord );
		for( int row_i = 0; row_i < rowField->mesh->el->nPoints; row_i++ ) {
			for( int col_j = 0; col_j < colField->mesh->el->nPoints; col_j++ ) {
				M[row_i*(2*colField->mesh->el->nNodes)+(2*col_j+0)] -= detJac*weight*constant*Ni[row_i]*GNix[0][col_j];
				M[row_i*(2*colField->mesh->el->nNodes)+(2*col_j+1)] -= detJac*weight*constant*Ni[row_i]*GNix[1][col_j];
			}
		}
	}
#endif
/**/
/* ...so does this but converges verrrrry slooooooowwwwwwwlllllllllllllyyyyyyyyyyyyyyyyy.... 
	double *coord, detJac, weight, **GNix, *Ni;
	for( int pt_i = 0; pt_i < rowField->mesh->el->nPoints; pt_i++ ) {
		coord  = rowField->mesh->el->quadPts[pt_i]->coord;
		weight = rowField->mesh->el->quadPts[pt_i]->weight;
		GNix   = rowField->mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		Ni     = colField->mesh->el->ShapeFuncsAtCoord( coord );
		for( int row_i = 0; row_i < rowField->mesh->el->nNodes; row_i++ ) {
			for( int col_j = 0; col_j < colField->mesh->el->nNodes; col_j++ ) {
				M[row_i*(2*colField->mesh->el->nNodes)+(2*col_j+0)] -= detJac*weight*constant*GNix[0][row_i]*Ni[col_j];
				M[row_i*(2*colField->mesh->el->nNodes)+(2*col_j+1)] -= detJac*weight*constant*GNix[1][row_i]*Ni[col_j];
			}
		}
	}
*/
	int row, col;
	double detJac, weight, *coord;
	Mesh* rowMesh = colField->mesh; /* transpose of the gradient matrix */
	Mesh* colMesh = rowField->mesh;
	Element* rowEl = rowMesh->el;
	Element* colEl = colMesh->el;
	int nRowNodes = rowEl->nNodes;
	int nColNodes = colEl->nNodes;
	int nRowDofs = colField->nDofs;
	int nColDofs = rowField->nDofs;
	int nRowEntries = nRowNodes*nRowDofs;
	double** GNix_row;
	double* Ni_col;
	double w;

	for( int pt_i = 0; pt_i < rowEl->nPoints; pt_i++ ) {
		coord  = rowEl->quadPts[pt_i]->coord;
		weight = rowEl->quadPts[pt_i]->weight;
		GNix_row = rowMesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		Ni_col = colEl->ShapeFuncsAtCoord( coord );
		w = detJac*weight*constant;
		for( int row_i = 0; row_i < nRowNodes; row_i++ ) {
			for( int dof_i = 0; dof_i < nRowDofs; dof_i++ ) {
				for( int col_j = 0; col_j < nColNodes; col_j++ ) {
					for( int dof_j = 0; dof_j < nColDofs; dof_j++ ) {
						row = row_i*nRowDofs + dof_i;
						col = col_j*nColDofs + dof_j;
						M[col*nRowEntries + row] += w*GNix_row[dof_i][row_i]*Ni_col[col_j];
					}
				}
			}
		}
	}
}

StokesPC::StokesPC( string _name, Field* _rowField, Field* _colField, double _constant, Field* _velocity, Field* _pressure ) : Operator( _name, _rowField, _colField, _constant ) {
	int nVelNodes, nVelDofs, nPresNodes;

	velocity = _velocity;
	pressure = _pressure;

	nVelNodes = velocity->mesh->el->nNodes;
	nVelDofs = velocity->nDofs;
	nPresNodes = pressure->mesh->el->nNodes;

	K = new double[nVelNodes*nVelDofs*nVelNodes*nVelDofs];
	G = new double[nVelNodes*nVelDofs*nPresNodes];
	KinvG = new double[nVelNodes*nVelDofs*nPresNodes];

	laplacian = new Laplacian( "stokesPC_laplacian", velocity, velocity, constant );
	gradient = new Gradient( "stokesPC_gradient", velocity, pressure, 1.0 );
}

StokesPC::~StokesPC() {
	delete[] K;
	delete[] G;
	delete[] KinvG;
	delete laplacian;
	delete gradient;
}

void StokesPC::AssembleElement( int el_i, double* M ) {
	int nVelNodesAndDofs = velocity->nDofs*velocity->mesh->el->nNodes;
	int nPresNodes = pressure->mesh->el->nNodes;

	memset( K, 0, nVelNodesAndDofs*nVelNodesAndDofs*sizeof(double) );
	memset( G, 0, nVelNodesAndDofs*nPresNodes*sizeof(double) );
	laplacian->AssembleElement( el_i, K );
	gradient->AssembleElement( el_i, G );

	for( int row_i = 0; row_i < nVelNodesAndDofs; row_i++ ) {
		for( int col_j = 0; col_j < nPresNodes; col_j++ ) {
			KinvG[row_i*nPresNodes + col_j] = G[row_i*nPresNodes + col_j] / K[row_i*nVelNodesAndDofs + row_i];
		}
	}

	for( int row_i = 0; row_i < nPresNodes; row_i++ ) {
		for( int col_j = 0; col_j < nPresNodes; col_j++ ) {
			for( int vel_i = 0; vel_i < nVelNodesAndDofs; vel_i++ ) {
				M[row_i*nPresNodes + col_j] += constant*G[vel_i*nPresNodes + col_j]*KinvG[vel_i*nPresNodes + row_i];
			}
		}
	}
}

BetaMatrix::BetaMatrix( string _name, Field* _rowField, Field* _colField, double _constant, double _f0, double _beta ) : Operator( _name, _rowField, _colField, _constant ) {
	f0 = _f0;
	beta = _beta;
}

BetaMatrix::~BetaMatrix() {}

void BetaMatrix::AssembleElement( int el_i, double* M ) {
	int row, col;
	double *coord, weight, gCoord[3], detJac, f, w, *Ni;
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
		w = detJac*weight*f*constant;
		for( int row_i = 0; row_i < nNodes; row_i++ ) {
			for( int col_j = 0; col_j < nNodes; col_j++ ) {
				row = row_i*2;
				col = col_j*2;
				M[(row+0)*nEntries + (col+1)] -= w*Ni[row_i]*Ni[col_j];
				M[(row+1)*nEntries + (col+0)] += w*Ni[row_i]*Ni[col_j];
			}
		}
	}
}

ConvectionMatrix::ConvectionMatrix( string _name, Field* _rowField, Field* _colField, double _constant, Field* _field ) : Operator( _name, _rowField, _colField, _constant ) {
        field = _field;
}

ConvectionMatrix::~ConvectionMatrix() {}

void ConvectionMatrix::AssembleElement( int el_i, double* M ) {
        int row, col;
        double weight, detJac, *Ni, **GNix, vel[2];
        Mesh* mesh = rowField->mesh;
        Element* el = mesh->el;
        int nNodes = el->nNodes;
        int nEntries = 2*nNodes;
	int* elNodes = mesh->ElNodes( el_i );

        for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
                weight = el->quadPts[pt_i]->weight;
                Ni = el->ShapeFuncs( pt_i );
                GNix = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		vel[0] = field->vals[elNodes[pt_i]][0];
		vel[1] = field->vals[elNodes[pt_i]][1];
                for( int row_i = 0; row_i < nNodes; row_i++ ) {
                        for( int col_j = 0; col_j < nNodes; col_j++ ) {
                                row = row_i*2;
                                col = col_j*2;
                                M[(row+0)*nEntries + (col+0)] -= detJac*weight*Ni[row_i]*(vel[0]*GNix[0][col_j] + vel[1]*GNix[1][col_j]);
                                M[(row+1)*nEntries + (col+1)] -= detJac*weight*Ni[row_i]*(vel[0]*GNix[0][col_j] + vel[1]*GNix[1][col_j]);
                        }
                }
        }
}
