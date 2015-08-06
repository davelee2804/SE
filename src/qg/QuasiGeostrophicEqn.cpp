#include <iostream>
#include <fstream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "LinSys.h"
#include "QuasiGeostrophicEqn.h"

using namespace std;
using std::string;

/*
Solve the QG eqn:
	D/Dt{zeta} = 0
as
	1: advect zeta 		=> zeta = omega - F.phi + beta.y
	2: solve for phi	=> (nabla^2 - F).phi = zeta - beta.y
	3: find omega 		=> omega = zeta - beta.y + F.phi
with linear shear
*/

QuasiGeostrophicEqn::QuasiGeostrophicEqn( Field* _phi, Field* _omega, double _F, double _beta, double _nu, double _r, double _dt, Field* _topo ) {
	DisRHS*		disRHS;
	RHSOp**		rhsOps;

	omega 	= _omega;
	phi	= _phi;

	F	= _F;
	beta 	= _beta;
	nu	= _nu;
	r	= _r;
	dt	= _dt;
	topo	= _topo;

	alpha		= new Field( "alpha", phi->mesh, 1, phi->bcs );
	alphaHat	= new Field( "alpha-hat", phi->mesh, 1, phi->bcs );
	phiHat		= new Field( "phi-hat", phi->mesh, 1, phi->bcs );
	omegaHat	= new Field( "omega-hat", phi->mesh, 1, phi->bcs );
	helmRhs1	= new Field( "helm-rhs-1", phi->mesh, 1, NULL );
	helmRhs2	= new Field( "helm-rhs-2", phi->mesh, 1, NULL );
	he1 		= new HelmholtzEqn( phiHat, helmRhs1, -F, -1.0, 1.0, 0 );
	he2 		= new HelmholtzEqn( phi, helmRhs2, -F, -1.0, 1.0, 0 );

	he1->Assemble();
	he2->Assemble();

	/* setup the dissipative terms solver */
	secondOrderAssembled = false;
	AssembleQMat( 1.0 );

	VecDuplicate( b, &f );

	disRHS = new DisRHS( "dissipation-rhs", phi->mesh, dt, alphaHat, omegaHat, r, nu );
	rhsOps = new RHSOp*[1];
	rhsOps[0] = disRHS;
	qRHS = new Vector( "q-rhs", alphaHat, f, rhsOps, 1 );
}

QuasiGeostrophicEqn::~QuasiGeostrophicEqn() {
	delete alpha;
	delete alphaHat;
	delete phiHat;
	delete omegaHat;
	delete helmRhs1;
	delete helmRhs2;
	delete he1;
	delete he2;

	delete qRHS;

	VecDestroy( b );
	VecDestroy( f );
	MatDestroy( M );
}

void QuasiGeostrophicEqn::Solve( Field* phiPrev, Field* omegaPrev ) {
	Field*	phiTemp 	= NULL;
	Field*	omegaTemp	= NULL;

	if( omegaPrev != NULL && phiPrev != NULL ) {
		if( !secondOrderAssembled ) {
			/* copy in the bcs */
			for( int node_i = 0; node_i < phi->mesh->nVertsTotal; node_i++ ) {
				phiHat->vals[node_i][0] = 1.5*phi->vals[node_i][0];
				omegaHat->vals[node_i][0] = 1.5*omega->vals[node_i][0];
			}
			MatDestroy( M );
			VecDestroy( b );
			AssembleQMat( 1.5 );
			secondOrderAssembled = true;
		}

		phiTemp 	= new Field( "phi-temp", phi->mesh, 1, NULL );
		omegaTemp	= new Field( "omega-temp", phi->mesh, 1, NULL );
		phiTemp->Copy( phi );
		omegaTemp->Copy( omega );

		SecondOrder( phiPrev, omegaPrev );

		phiPrev->Copy( phiTemp );
		omegaPrev->Copy( omegaTemp );
		delete phiTemp;
		delete omegaTemp;
	}
	else {
		FirstOrder();
	}
}

void QuasiGeostrophicEqn::FirstOrder() {
	Field*		vel	= GenVel( phi );
	Advector*	adv	= new Advector( alpha, vel );

	GenPV( alpha, phi, omega );

	adv->Advect( dt );
	for( int vert_i = 0; vert_i < alpha->mesh->nVertsTotal; vert_i++ ) {
		helmRhs1->vals[vert_i][0] = adv->fieldSL->vals[vert_i][0] - beta*alpha->mesh->verts[vert_i][1];
		if( topo ) {
			helmRhs1->vals[vert_i][0] -= topo->vals[vert_i][0];
		}
	}
	he1->Solve( "helm_" );
	for( int vert_i = 0; vert_i < alpha->mesh->nVertsTotal; vert_i++ ) {
		omegaHat->vals[vert_i][0] = adv->fieldSL->vals[vert_i][0] + F*phiHat->vals[vert_i][0] - beta*alpha->mesh->verts[vert_i][1];
		if( topo ) {
			omegaHat->vals[vert_i][0] -= topo->vals[vert_i][0];
		}
	}

	alphaHat->Copy( adv->fieldSL );
	SolveQ();

	for( int vert_i = 0; vert_i < alpha->mesh->nVertsTotal; vert_i++ ) {
		helmRhs2->vals[vert_i][0] = alpha->vals[vert_i][0] - beta*alpha->mesh->verts[vert_i][1];
		if( topo ) {
			helmRhs2->vals[vert_i][0] -= topo->vals[vert_i][0];
		}
	}
	he2->Solve( "helm_" );
	for( int vert_i = 0; vert_i < alpha->mesh->nVertsTotal; vert_i++ ) {
		omega->vals[vert_i][0] = alpha->vals[vert_i][0] + F*phi->vals[vert_i][0] - beta*alpha->mesh->verts[vert_i][1];
		if( topo ) {
			omega->vals[vert_i][0] -= topo->vals[vert_i][0];
		}
	}

	delete vel;
	delete adv;
}

void QuasiGeostrophicEqn::SecondOrder( Field* phiPrev, Field* omegaPrev ) {
	Field*		vel		= GenVel( phi );
	Field*		velPrev		= GenVel( phiPrev );
	Field*		alphaPrev	= new Field( "alpha-prev", phi->mesh, 1, NULL );
	Advector*	adv		= new Advector( alpha, vel, alphaPrev, velPrev );
	double		a		= 1.5;

	GenPV( alpha, phi, omega );
	GenPV( alphaPrev, phiPrev, omegaPrev );

	adv->Advect( dt );
	for( int vert_i = 0; vert_i < alpha->mesh->nVertsTotal; vert_i++ ) {
		helmRhs1->vals[vert_i][0] = adv->fieldSL->vals[vert_i][0] - a*beta*alpha->mesh->verts[vert_i][1];
		if( topo ) {
			helmRhs1->vals[vert_i][0] -= a*topo->vals[vert_i][0];
		}
	}
	he1->Solve( "helm_" );
	for( int vert_i = 0; vert_i < alpha->mesh->nVertsTotal; vert_i++ ) {
		omegaHat->vals[vert_i][0] = adv->fieldSL->vals[vert_i][0] + F*phiHat->vals[vert_i][0] - a*beta*alpha->mesh->verts[vert_i][1];
		if( topo ) {
			omegaHat->vals[vert_i][0] -= a*topo->vals[vert_i][0];
		}
	}

	alphaHat->Copy( adv->fieldSL );
	SolveQ();

	for( int vert_i = 0; vert_i < alpha->mesh->nVertsTotal; vert_i++ ) {
		helmRhs2->vals[vert_i][0] = alpha->vals[vert_i][0] - beta*alpha->mesh->verts[vert_i][1];
		if( topo ) {
			helmRhs2->vals[vert_i][0] -= topo->vals[vert_i][0];
		}
	}
	he2->Solve( "helm_" );
	for( int vert_i = 0; vert_i < alpha->mesh->nVertsTotal; vert_i++ ) {
		omega->vals[vert_i][0] = alpha->vals[vert_i][0] + F*phi->vals[vert_i][0] - beta*alpha->mesh->verts[vert_i][1];
		if( topo ) {
			omega->vals[vert_i][0] -= topo->vals[vert_i][0];
		}
	}

	delete vel;
	delete velPrev;
	delete alphaPrev;
	delete adv;
}

Field* QuasiGeostrophicEqn::GenVel( Field* psi ) {
	Field*		vel 	= new Field( "vel", psi->mesh, 2, NULL );
	double**	dp	= new double*[1];

	dp[0] = new double[2];

	for( int vert_i = 0; vert_i < psi->mesh->nVertsTotal; vert_i++ ) {
		psi->InterpDerivsGlobal( phi->mesh->verts[vert_i], dp );
		vel->vals[vert_i][0] = -1.0*dp[0][1];
		vel->vals[vert_i][1] = +1.0*dp[0][0];
	}

	delete[] dp[0];
	delete[] dp;

	return vel;
}

void QuasiGeostrophicEqn::GenPV( Field* pv, Field* psi, Field* vort ) {
	for( int vert_i = 0; vert_i < psi->mesh->nVertsTotal; vert_i++ ) {
		pv->vals[vert_i][0] = vort->vals[vert_i][0] - F*psi->vals[vert_i][0] + beta*psi->mesh->verts[vert_i][1];
		if( topo ) {
			pv->vals[vert_i][0] += topo->vals[vert_i][0];
		}
	}
}

void QuasiGeostrophicEqn::AssembleQMat( double a ) {
	int		size	= alpha->mesh->nVertsTotal - alpha->bcs->size[0];
	Vec		x;
	Vector*		sol;
	Vector*		rhs;
	Matrix*		mat;
	Operator**	ops;
	MassMatrix*	massMat;

	MatCreate( MPI_COMM_WORLD, &M );
	MatSetSizes( M, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( M, MATSEQAIJ );
	MatSetFromOptions( M );
	MatSeqAIJSetPreallocation( M, 4*alpha->mesh->el->nNodes, PETSC_NULL );
	MatZeroEntries( M );

	VecCreate( MPI_COMM_WORLD, &b );
	VecSetSizes( b, size, PETSC_DETERMINE );
	VecSetType( b, VECSEQ );
	VecSetFromOptions( b );
	VecSetOption( b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
	VecZeroEntries( b );

	VecDuplicate( b, &x );
	VecZeroEntries( x );

	massMat = new MassMatrix( "mass-mat", alpha, alpha, a );
	ops = new Operator*[1];
	ops[0] = massMat;

	sol = new Vector( "sol", alpha, x, NULL, 0 );
	rhs = new Vector( "rhs", alpha, b, NULL, 0 );
	mat = new Matrix( "mat", M, sol, sol, rhs, ops, 1 );
	mat->Assemble();
	MatAssemblyBegin( M, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( M, MAT_FINAL_ASSEMBLY );

	VecDestroy( x );
	delete sol;
	delete rhs;
	delete mat;
}

void QuasiGeostrophicEqn::SolveQ() {
	Vector* sol;
	Vec	x;
	KSP	ksp;

	VecDuplicate( f, &x );
	VecZeroEntries( x );
	VecZeroEntries( f );
	sol = new Vector( "sol", alpha, x, NULL, 0 );
	qRHS->Assemble();
	VecAXPY( f, 1.0, b );
	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, M, M, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "q_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	sol->UpdateField();

	KSPDestroy( ksp );
	VecDestroy( x );
	delete sol;
}

/* RHS operator */
DisRHS::DisRHS( string _name, Mesh* _mesh, double _constant, Field* _field, Field* _omega, double _r, double _nu ) : RHSOp( _name, _mesh, _constant, _field ) {
	omega	= _omega;
	r	= _r;
	nu	= _nu;

	gW = new double*[1];
	gW[0] = new double[2];
}

DisRHS::~DisRHS() {
	delete[] gW[0];
	delete[] gW;
}

void DisRHS::AssembleElement( int el_i, double* rhs ) {
	double	*Ni, **GNx, *coord, detJac, weight, c, q, w, gCoord[2];

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord	= mesh->el->quadPts[pt_i]->coord;
		weight	= mesh->el->quadPts[pt_i]->weight;
		Ni 	= mesh->el->ShapeFuncs( pt_i );
		GNx	= mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );

		c = detJac*weight*constant;

		mesh->LocalToGlobal( coord, el_i, gCoord );
		field->InterpLocal( el_i, coord, &q );
		omega->InterpLocal( el_i, coord, &w );
		omega->InterpDerivsGlobal( gCoord, gW );

		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] += detJac*weight*q*Ni[node_i];
			rhs[node_i] -= c*r*w*Ni[node_i];
			rhs[node_i] -= c*nu*( GNx[0][node_i]*gW[0][0] + GNx[1][node_i]*gW[0][1] );
		}
	}
}
