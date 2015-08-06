#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Legendre.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "RHSOp.h"

using namespace std;
using std::string;

RHSOp::RHSOp( string _name, Mesh* _mesh, double _constant, Field* _field ) {
	name     = _name;
	mesh     = _mesh;
	constant = _constant;
	field    = _field;
}

RHSOp::~RHSOp() {}

FieldRHS::FieldRHS( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {}

FieldRHS::~FieldRHS() {}

void FieldRHS::AssembleElement( int el_i, double* rhs ) {
	double* Ni;
	double* coord, detJac, weight;
	double f[3];
	double w;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		Ni = mesh->el->ShapeFuncs( pt_i );
		detJac = mesh->DetJac( el_i, pt_i );
		field->InterpLocal( el_i, coord, f );
		w = detJac*weight*constant;
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			for( int dof_i = 0; dof_i < field->nDofs; dof_i++ ) {
				rhs[node_i*field->nDofs + dof_i] += w*f[dof_i]*Ni[node_i];
			}
		}
	}
}

GradFieldRHS::GradFieldRHS( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {}

GradFieldRHS::~GradFieldRHS() {}

void GradFieldRHS::AssembleElement( int el_i, double* rhs ) {
	double *coord, detJac, weight, **GNix, f;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpLocal( el_i, coord, &f );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[2*node_i+0] -= detJac*weight*constant*GNix[0][node_i]*f;
			rhs[2*node_i+1] -= detJac*weight*constant*GNix[1][node_i]*f;
		}
	}
}

DivVelRHS::DivVelRHS( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {
	if( field->nDofs != 2 ) { cerr << "div-vel-rhs-op field is not a vector field!\n"; exit(0); }
}

DivVelRHS::~DivVelRHS() {}

void DivVelRHS::AssembleElement( int el_i, double* rhs ) {
	double *coord, detJac, weight, u[2], **GNix;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpLocal( el_i, coord, u );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] -= detJac*weight*constant*( u[0]*GNix[0][node_i] + u[1]*GNix[1][node_i] );
		}
	}
}

GradDotFieldRHS::GradDotFieldRHS( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {}

GradDotFieldRHS::~GradDotFieldRHS() {}

void GradDotFieldRHS::AssembleElement( int el_i, double* rhs ) {
	double *coord, detJac, weight, f[2], **GNix;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpLocal( el_i, coord, f );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] += detJac*weight*( f[0]*GNix[0][node_i] + f[1]*GNix[1][node_i] );
		}
	}
}

ConvectionRHS::ConvectionRHS( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {
        du = new double*[2]; du[0] = new double[2]; du[1] = new double[2];
}

ConvectionRHS::~ConvectionRHS() {
        delete[] du[0]; delete[] du[1]; delete[] du;
}

void ConvectionRHS::AssembleElement( int el_i, double* rhs ) {
        double *coord, weight, detJac, *Ni, u[2], gCoord[2];

        for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
                coord  = mesh->el->quadPts[pt_i]->coord;
                weight = mesh->el->quadPts[pt_i]->weight;
                detJac = mesh->DetJac( el_i, pt_i );
                Ni     = mesh->el->ShapeFuncs( pt_i );
                mesh->LocalToGlobal( coord, el_i, gCoord );
                field->InterpLocal( el_i, coord, u );
                field->InterpDerivsGlobal( gCoord, du );
                for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
                        rhs[2*node_i+0] += detJac*weight*constant*Ni[node_i]*( u[0]*du[0][0] + u[1]*du[0][1] );
                        rhs[2*node_i+1] += detJac*weight*constant*Ni[node_i]*( u[0]*du[1][0] + u[1]*du[1][1] );
                }
        }
}

StressTensorRHS::StressTensorRHS( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {
        du = new double*[2]; du[0] = new double[2]; du[1] = new double[2];
}

StressTensorRHS::~StressTensorRHS() {
        delete[] du[0]; delete[] du[1]; delete[] du;
}

void StressTensorRHS::AssembleElement( int el_i, double* rhs ) {
	double weight, detJac, **GNix;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[2*node_i+0] -= detJac*weight*constant*( 2.0*GNix[0][node_i]*du[0][0] + GNix[1][node_i]*( du[0][1] + du[1][0] ) );
			rhs[2*node_i+1] -= detJac*weight*constant*( 2.0*GNix[1][node_i]*du[1][1] + GNix[0][node_i]*( du[0][1] + du[1][0] ) );
		}
	}
}

DivHeightVelRHS::DivHeightVelRHS( string _name, Mesh* _mesh, double _constant, Field* _field, Field* _velocity ) : 
RHSOp( _name, _mesh, _constant, _field ) 
{
	velocity = _velocity;
}

DivHeightVelRHS::~DivHeightVelRHS() {}

void DivHeightVelRHS::AssembleElement( int el_i, double* rhs ) {
	double *coord, detJac, weight, **GNix, h, v[2], a;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpLocal( el_i, coord, &h );
		velocity->InterpLocal( el_i, coord, v );
		a = detJac*weight*constant*h;
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] -= a*( GNix[0][node_i]*v[0] + GNix[1][node_i]*v[1] );
		}
	}
}

LaplacianRHS::LaplacianRHS( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {}

LaplacianRHS::~LaplacianRHS() {}

void LaplacianRHS::AssembleElement( int el_i, double* rhs ) {
	double *coord, detJac, weight, **GNix, **df, gCoord[2];
	df = new double*[1];
	df[0] = new double[2];

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		mesh->LocalToGlobal( coord, el_i, gCoord );
		field->InterpDerivsGlobal( gCoord, df );
		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] -= detJac*weight*constant*GNix[0][node_i]*df[0][0];
		}
	}

	delete[] df[0];
	delete[] df;
}

BiharmonicSurfaceRHS::BiharmonicSurfaceRHS( string _name, Mesh* _mesh, double _constant, Field* _field ) : RHSOp( _name, _mesh, _constant, _field ) {
	g2 = new double*[1];
	g2[0] = new double[3];
}

BiharmonicSurfaceRHS::~BiharmonicSurfaceRHS() {
	delete[] g2[0];
	delete[] g2;
}

void BiharmonicSurfaceRHS::AssembleElement( int el_i, double* rhs ) {
	Legendre*	el 	= (Legendre*)mesh->el;
	double 		norm, detJac, weight, **GNix, *coord, gCoord[2];
	int		N	= el->N;

	for( int pt_i = 0; pt_i < el->nPoints; pt_i++ ) {
		coord = el->quadPts[pt_i]->coord;
		mesh->LocalToGlobal( coord, el_i, gCoord );
		GNix = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );
		field->InterpSecondDerivs( el_i, pt_i, g2 );
		/* bottom */
		if( abs(gCoord[1] - mesh->min[1]) < 1.0e-6 ) {
			norm   = -1.0;
			detJac = 0.5*mesh->dx[1];
			weight = el->weights[pt_i%(el->N+1)];
			/* 2nd derivative surface integral */
			for( int node_i = 0; node_i <= N; node_i++ ) {
				rhs[node_i] -= detJac*weight*norm*GNix[1][node_i]*(g2[0][0]+g2[0][1]);
			}
		}
		/* top */
		if( abs(gCoord[1] - mesh->max[1]) < 1.0e-6 ) {
			norm   = +1.0;
			detJac = 0.5*mesh->dx[1];
			weight = el->weights[pt_i%(el->N+1)];
			/* 2nd derivative surface integral */
			for( int node_i = 0; node_i <= N; node_i++ ) {
				rhs[(N+1)*N+node_i] -= detJac*weight*norm*GNix[1][(N+1)*N+node_i]*(g2[0][0]+g2[0][1]);
			}
		}
		/* left */
		if( abs(gCoord[0] - mesh->min[0]) < 1.0e-6 ) {
			norm   = -1.0;
			detJac = 0.5*mesh->dx[0];
			weight = el->weights[pt_i/(el->N+1)];
			/* 2nd derivative surface integral */
			for( int node_i = 0; node_i <= N; node_i++ ) {
				rhs[node_i*(N+1)] -= detJac*weight*norm*GNix[0][node_i*(N+1)]*(g2[0][0]+g2[0][1]);
			}
		}
		/* right */
		if( abs(gCoord[0] - mesh->max[0]) < 1.0e-6 ) {
			norm   = +1.0;
			detJac = 0.5*mesh->dx[0];
			weight = el->weights[pt_i/(el->N+1)];
			/* 2nd derivative surface integral */
			for( int node_i = 0; node_i <= N; node_i++ ) {
				rhs[node_i*(N+1)+N] -= detJac*weight*norm*GNix[0][node_i*(N+1)+N]*(g2[0][0]+g2[0][1]);
			}
		}
	}
}
