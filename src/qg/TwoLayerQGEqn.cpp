#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "LinSys.h"
#include "TwoLayerQGEqn.h"

using namespace std;
using std::string;

TwoLayerQGEqn::TwoLayerQGEqn( Field* _psi1, Field* _omega1, Field* _psi2, Field* _omega2, double _F1, double _F2, double _beta, double _E2, double _r, double _nu, Field* tau_e, double _dt ) {
	DissipationRHS*	d1RHS;
	DissipationRHS*	d2RHS;
	RHSOp**		rhsOp1;
	RHSOp**		rhsOp2;

	omega1	= _omega1;
	omega2	= _omega2;
	psi1	= _psi1;
	psi2	= _psi2;
	F1	= _F1;
	F2	= _F2;
	beta	= _beta;
	E2	= _E2;
	r	= _r;
	nu	= _nu;
	dt	= _dt;

	tau1 	= new Field( "tau-1", psi1->mesh, 1, psi1->bcs );
	chi1 	= new Field( "chi-1", psi1->mesh, 1, psi1->bcs );
	tau2 	= new Field( "tau-2", psi1->mesh, 1, psi1->bcs );
	chi2	= new Field( "chi-2", psi1->mesh, 1, psi1->bcs );

	psi1Hat = new Field( "psi-1-hat", psi1->mesh, 1, psi1->bcs );
	psi2Hat = new Field( "psi-1-hat", psi1->mesh, 1, psi1->bcs );
	omega1Hat = new Field( "omega-1-hat", psi1->mesh, 1, psi1->bcs );
	omega2Hat = new Field( "omega-1-hat", psi1->mesh, 1, psi1->bcs );

	/* set the bcs */
	for( int vert_i = 0; vert_i < psi1->mesh->nVertsTotal; vert_i++ ) {
		if( psi1->bcs->IsBCNode( vert_i ) ) {
			tau1->vals[vert_i][0] = psi1->vals[vert_i][0] - psi2->vals[vert_i][0];
			chi1->vals[vert_i][0] = psi1->vals[vert_i][0] + psi2->vals[vert_i][0];
			tau2->vals[vert_i][0] = psi1->vals[vert_i][0] - psi2->vals[vert_i][0];
			chi2->vals[vert_i][0] = psi1->vals[vert_i][0] + psi2->vals[vert_i][0];
		}
	}

	helmRhs1 = new Field( "helm-rhs-1", psi1->mesh, 1, NULL );
	fishRhs1 = new Field( "fish-rhs-1", psi1->mesh, 1, NULL );
	helmRhs2 = new Field( "helm-rhs-2", psi1->mesh, 1, NULL );
	fishRhs2 = new Field( "fish-rhs-2", psi1->mesh, 1, NULL );

	helmEqn1 = new HelmholtzEqn( tau1, helmRhs1, -(F1 + F2), -1.0, 1.0, 0 );
	helmEqn2 = new HelmholtzEqn( tau2, helmRhs2, -(F1 + F2), -1.0, 1.0, 0 );
	fishEqn1 = new PoissonEqn( chi1, fishRhs1, -1.0 );
	fishEqn2 = new PoissonEqn( chi2, fishRhs2, -1.0 );

	/* generate the potential vorticity mass matrices */
	secondOrderAssembled = false;
	alpha1 = new Field( "q1", psi1->mesh, 1, psi1->bcs );
	alpha2 = new Field( "q2", psi2->mesh, 1, psi2->bcs );
	alpha1Hat = new Field( "q1-hat", psi1->mesh, 1, psi1->bcs );
	alpha2Hat = new Field( "q2-hat", psi2->mesh, 1, psi2->bcs );
	GenPV( alpha1, omega1, psi1, psi2, F1 );
	GenPV( alpha2, omega2, psi2, psi1, F2 );
	AssembleQMat( &M1, &b1, alpha1, 1.0 );
	AssembleQMat( &M2, &b2, alpha2, 1.0 );
	VecDuplicate( b1, &f1 );
	VecDuplicate( b2, &f2 );

	d1RHS = new DissipationRHS( "d1", psi1->mesh, dt, alpha1Hat, tau1, omega1Hat, +r*F1, 0.0, tau_e, nu );
	rhsOp1 = new RHSOp*[1];
	rhsOp1[0] = d1RHS;
	q1RHS = new Vector( "q1-rhs", alpha1, f1, rhsOp1, 1 );

	d2RHS = new DissipationRHS( "d2", psi1->mesh, dt, alpha2Hat, tau1, omega2Hat, -r*F2, E2,  tau_e, nu );
	rhsOp2 = new RHSOp*[1];
	rhsOp2[0] = d2RHS;
	q2RHS = new Vector( "q2-rhs", alpha2, f2, rhsOp2, 1 );
}

TwoLayerQGEqn::~TwoLayerQGEqn() {
	delete tau1;
	delete chi1;
	delete tau2;
	delete chi2;

	delete helmRhs1;
	delete fishRhs1;
	delete helmRhs2;
	delete fishRhs2;

	delete fishEqn1;
	delete helmEqn1;
	delete fishEqn2;
	delete helmEqn2;

	delete phi1Hat;
	delete phi2Hat;
	delete omega1Hat;
	delete omega2Hat;

	delete alpha1;
	delete alpha2;
	delete alpha1Hat;
	delete alpha2Hat;

	MatDestroy( M1 );
	MatDestroy( M2 );
	VecDestroy( b1 );
	VecDestroy( b2 );
	VecDestroy( f1 );
	VecDestroy( f2 );
}

void TwoLayerQGEqn::Solve( Field* psi1Prev, Field* omega1Prev, Field* psi2Prev, Field* omega2Prev ) {
	Field*		psi1Temp 	= NULL;
	Field*		psi2Temp	= NULL;
	Field*		omega1Temp	= NULL;
	Field*		omega2Temp	= NULL;

	if( psi1Prev != NULL && omega1Prev != NULL && psi2Prev != NULL && omega2Prev != NULL ) {
		if( !secondOrderAssembled ) {
			/* copy in the bcs */
			for( int node_i = 0; node_i < psi1->mesh->nVertsTotal; node_i++ ) {
				psi1Hat->vals[node_i][0] = 1.5*psi1->vals[node_i][0];
				psi2Hat->vals[node_i][0] = 1.5*psi2->vals[node_i][0];
				omega1Hat->vals[node_i][0] = 1.5*omega1->vals[node_i][0];
				omega2Hat->vals[node_i][0] = 1.5*omega2->vals[node_i][0];
			}
			MatDestroy( M1 );
			MatDestroy( M2 );
			VecDestroy( b1 );
			VecDestroy( b2 );
			AssembleQMat( &M1, &b1, alpha1, 1.5 );
			AssembleQMat( &M2, &b2, alpha2, 1.5 );
			secondOrderAssembled = true;
		}

		psi1Temp 	= new Field( "phi-1-temp", psi1->mesh, 1, NULL );
		psi2Temp 	= new Field( "phi-2-temp", psi2->mesh, 1, NULL );
		omega1Temp	= new Field( "omega-1-temp", omega1->mesh, 1, NULL );
		omega2Temp	= new Field( "omega-2-temp", omega2->mesh, 1, NULL );

		psi1Temp->Copy( psi1 );
		psi2Temp->Copy( psi2 );
		omega1Temp->Copy( omega1 );
		omega2Temp->Copy( omega2 );

		SecondOrder( psi1Prev, omega1Prev, psi2Prev, omega2Prev );

		psi1Prev->Copy( psi1Temp );
		psi2Prev->Copy( psi2Temp );
		omega1Prev->Copy( omega1Temp );
		omega2Prev->Copy( omega2Temp );

		delete psi1Temp;
		delete psi2Temp;
		delete omega1Temp;
		delete omega2Temp;
	}
	else {
		FirstOrder();
	}
}

void TwoLayerQGEqn::FirstOrder() {
	Field*		vel1		= GenVel( psi1, "vel-1" );
	Field*		vel2		= GenVel( psi2, "vel-2" );
	Advector*	adv1		= new Advector( alpha1, vel1 );
	Advector*	adv2		= new Advector( alpha2, vel2 );
	double		y;

	/* 1. advection of potential vorticity */
	GenPV( alpha1, omega1, psi1, psi2, F1 );
	GenPV( alpha2, omega2, psi2, psi1, F2 );
	
	adv1->Advect( dt );
	adv2->Advect( dt );
	alpha1Hat->Copy( adv1->fieldSL );
	alpha2Hat->Copy( adv2->fieldSL );

	/* 2. helmholtz solve for the baroclinic stream function */
	for( int node_i = 0; node_i < helmRhs1->mesh->nVertsTotal; node_i++ ) {
		helmRhs1->vals[node_i][0] = alpha1Hat->vals[node_i][0] - alpha2Hat->vals[node_i][0];
	}
	helmEqn1->Solve( "helm_" );

	/* 3. poisson solve for the barotropic stream function */
	for( int node_i = 0; node_i < fishRhs1->mesh->nVertsTotal; node_i++ ) {
		y = fishRhs1->mesh->verts[node_i][1];
		fishRhs1->vals[node_i][0] = alpha1Hat->vals[node_i][0] + alpha2Hat->vals[node_i][0] + 
					    (F1 - F2)*tau1->vals[node_i][0] - 2.0*beta*y;
	}
	fishEqn1->Solve( false );

	/* 4. update the stream functions and vorticities */
	for( int node_i = 0; node_i < psi1->mesh->nVertsTotal; node_i++ ) {
		y = psi1->mesh->verts[node_i][1];
		psi1Hat->vals[node_i][0] = 0.5*( tau1->vals[node_i][0] + chi1->vals[node_i][0] );
		psi2Hat->vals[node_i][0] = 0.5*( chi1->vals[node_i][0] - tau1->vals[node_i][0] );
		omega1Hat->vals[node_i][0] = alpha1Hat->vals[node_i][0] - F1*( psi2Hat->vals[node_i][0] - psi1Hat->vals[node_i][0] ) - beta*y;
		omega2Hat->vals[node_i][0] = alpha2Hat->vals[node_i][0] - F2*( psi1Hat->vals[node_i][0] - psi2Hat->vals[node_i][0] ) - beta*y;
	}

	/* 5. dissipative terms */
	SolveQ( M1, f1, b1, q1RHS, alpha1 );

	/* 6. helmholtz solve for the baroclinic stream function */
	for( int node_i = 0; node_i < helmRhs2->mesh->nVertsTotal; node_i++ ) {
		helmRhs2->vals[node_i][0] = alpha1->vals[node_i][0] - alpha2->vals[node_i][0];
	}
	helmEqn2->Solve( "helm_" );

	/* 7. poisson solve for the barotropic stream function */
	for( int node_i = 0; node_i < fishRhs2->mesh->nVertsTotal; node_i++ ) {
		y = fishRhs2->mesh->verts[node_i][1];
		fishRhs2->vals[node_i][0] = alpha1->vals[node_i][0] + alpha2->vals[node_i][0] + 
					    (F1 - F2)*tau2->vals[node_i][0] - 2.0*beta*y;
	}
	fishEqn2->Solve( false );

	/* 4. update the stream functions and vorticities */
	for( int node_i = 0; node_i < psi1->mesh->nVertsTotal; node_i++ ) {
		y = psi1->mesh->verts[node_i][1];
		psi1->vals[node_i][0] = 0.5*( tau2->vals[node_i][0] + chi2->vals[node_i][0] );
		psi2->vals[node_i][0] = 0.5*( chi2->vals[node_i][0] - tau2->vals[node_i][0] );
		omega1->vals[node_i][0] = alpha1->vals[node_i][0] - F1*( psi2->vals[node_i][0] - psi1->vals[node_i][0] ) - beta*y;
		omega2->vals[node_i][0] = alpha2->vals[node_i][0] - F2*( psi1->vals[node_i][0] - psi2->vals[node_i][0] ) - beta*y;
	}

	delete adv1;
	delete adv2;
	delete vel1;
	delete vel2;
}

void TwoLayerQGEqn::SecondOrder( Field* psi1Prev, Field* omega1Prev, Field* psi2Prev, Field* omega2Prev ) {
	Field*		vel1		= GenVel( psi1, "vel-1" );
	Field*		vel2		= GenVel( psi2, "vel-2" );
	Field*		vel1Prev 	= GenVel( psi1Prev, "vel-1-prev" );
	Field*		vel2Prev 	= GenVel( psi2Prev, "vel-2-prev" );
	Field*		alpha1Prev	= new Field( "alpha-1-prev", psi1->mesh, 1, NULL );
	Field*		alpha2Prev	= new Field( "alpha-2-prev", psi2->mesh, 1, NULL );
	Advector*	adv1		= new Advector( alpha1, vel1, alpha1Prev, vel1Prev );
	Advector*	adv2		= new Advector( alpha2, vel2, alpha2Prev, vel2Prev );
	double		a		= 1.5;
	double		y;

	/* 1. advection of potential vorticity */
	GenPV( alpha1, omega1, psi1, psi2, F1 );
	GenPV( alpha2, omega2, psi2, psi1, F2 );
	GenPV( alpha1Prev, omega1Prev, psi1Prev, psi2Prev, F1 );
	GenPV( alpha2Prev, omega2Prev, psi2Prev, psi1Prev, F2 );

	adv1->Advect( dt );
	adv2->Advect( dt );
	alpha1Hat->Copy( adv1->fieldSL );
	alpha2Hat->Copy( adv2->fieldSL );

	/* 2. helmholtz solve for the baroclinic stream function */
	for( int node_i = 0; node_i < helmRhs1->mesh->nVertsTotal; node_i++ ) {
		helmRhs1->vals[node_i][0] = alpha1Hat->vals[node_i][0] - alpha2Hat->vals[node_i][0];
	}
	helmEqn1->Solve( "helm_" );

	/* 3. poisson solve for the barotropic stream function */
	for( int node_i = 0; node_i < fishRhs1->mesh->nVertsTotal; node_i++ ) {
		y = fishRhs1->mesh->verts[node_i][1];
		fishRhs1->vals[node_i][0] = alpha1Hat->vals[node_i][0] + alpha2Hat->vals[node_i][0] + (F1 - F2)*tau1->vals[node_i][0] - a*2.0*beta*y;
	}
	fishEqn1->Solve( false );

	/* 4. update the stream functions and vorticities */
	for( int node_i = 0; node_i < psi1->mesh->nVertsTotal; node_i++ ) {
		y = psi1->mesh->verts[node_i][1];
		psi1Hat->vals[node_i][0] = 0.5*( tau1->vals[node_i][0] + chi1->vals[node_i][0] );
		psi2Hat->vals[node_i][0] = 0.5*( chi1->vals[node_i][0] - tau1->vals[node_i][0] );
		omega1Hat->vals[node_i][0] = alpha1Hat->vals[node_i][0] - F1*( psi2Hat->vals[node_i][0] - psi1Hat->vals[node_i][0] ) - a*beta*y;
		omega2Hat->vals[node_i][0] = alpha2Hat->vals[node_i][0] - F2*( psi1Hat->vals[node_i][0] - psi2Hat->vals[node_i][0] ) - a*beta*y;
	}

	/* 5. dissipative terms */
	SolveQ( M1, f1, b1, q1RHS, alpha1 );
	SolveQ( M2, f2, b1, q2RHS, alpha2 );

	/* 6. helmholtz solve for the baroclinic stream function */
	for( int node_i = 0; node_i < helmRhs2->mesh->nVertsTotal; node_i++ ) {
		helmRhs2->vals[node_i][0] = alpha1->vals[node_i][0] - alpha2->vals[node_i][0];
	}
	helmEqn2->Solve( "helm_" );

	/* 7. poisson solve for the barotropic stream function */
	for( int node_i = 0; node_i < fishRhs2->mesh->nVertsTotal; node_i++ ) {
		y = fishRhs2->mesh->verts[node_i][1];
		fishRhs2->vals[node_i][0] = alpha1->vals[node_i][0] + alpha2->vals[node_i][0] + (F1 - F2)*tau2->vals[node_i][0] - 2.0*beta*y;
	}
	fishEqn2->Solve( false );

	/* 4. update the stream functions and vorticities */
	for( int node_i = 0; node_i < psi1->mesh->nVertsTotal; node_i++ ) {
		y = psi1->mesh->verts[node_i][1];
		psi1->vals[node_i][0] = 0.5*( tau2->vals[node_i][0] + chi2->vals[node_i][0] );
		psi2->vals[node_i][0] = 0.5*( chi2->vals[node_i][0] - tau2->vals[node_i][0] );
		omega1->vals[node_i][0] = alpha1->vals[node_i][0] - F1*( psi2->vals[node_i][0] - psi1->vals[node_i][0] ) - beta*y;
		omega2->vals[node_i][0] = alpha2->vals[node_i][0] - F2*( psi1->vals[node_i][0] - psi2->vals[node_i][0] ) - beta*y;
	}

	delete adv1;
	delete adv2;
	delete vel1;
	delete vel2;
	delete vel1Prev;
	delete vel2Prev;
	delete alpha1Prev;
	delete alpha2Prev;
}

void TwoLayerQGEqn::GenPV( Field* alpha, Field* omega, Field* psi_i, Field* psi_j, double F ) {
	double	y;

	for( int node_i = 0; node_i < omega->mesh->nVertsTotal; node_i++ ) {
		y = omega->mesh->verts[node_i][1];
		alpha->vals[node_i][0] = omega->vals[node_i][0] + F*( psi_j->vals[node_i][0] - psi_i->vals[node_i][0] ) + beta*y;
	}
}

Field* TwoLayerQGEqn::GenVel( Field* psi, string name ) {
	double**	dp;
	Field*		vel 	= new Field( name, psi->mesh, 2, NULL );

	dp = new double*[1];
	dp[0] = new double[2];

	for( int node_i = 0; node_i < psi->mesh->nVertsTotal; node_i++ ) {
		psi->InterpDerivsGlobal( psi->mesh->verts[node_i], dp );
		vel->vals[node_i][0] = -dp[0][1];
		vel->vals[node_i][1] = +dp[0][0];
	}

	delete[] dp[0];
	delete[] dp;

	return vel;
}

void TwoLayerQGEqn::AssembleQMat( Mat* M, Vec* b, Field* qi, double a ) {
	int 		size 	= qi->mesh->nVertsTotal - qi->bcs->size[0];
	Vec		x;
	Vector*		sol;
	Vector*		rhs;
	Matrix*		mat;
	Operator**	ops	= new Operator*[1];
	MassMatrix*	massMat;

	MatCreate( MPI_COMM_WORLD, M );
	MatSetSizes( *M, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( *M, MATSEQAIJ );
	MatSetFromOptions( *M );
	MatSeqAIJSetPreallocation( *M, 4*qi->mesh->el->nNodes, PETSC_NULL );
	MatZeroEntries( *M );

	VecCreate( MPI_COMM_WORLD, &x );
	VecSetSizes( x, size, PETSC_DETERMINE );
	VecSetType( x, VECSEQ );
	VecSetFromOptions( x );
	VecSetOption( x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
	VecZeroEntries( x );

	VecDuplicate( x, b );
	VecZeroEntries( *b );

	massMat = new MassMatrix( "mass-mat", qi, qi, a );
	ops[0] = massMat;

	sol = new Vector( "sol", qi, x, NULL, 0 );
	rhs = new Vector( "rhs", qi, *b, NULL, 0 );
	mat = new Matrix( "q-mat", *M, sol, sol, rhs, ops, 1 );
	mat->Assemble();
	MatAssemblyBegin( *M, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( *M, MAT_FINAL_ASSEMBLY );

	VecDestroy( x );
	delete sol;
	delete rhs;
	delete mat;
}

void TwoLayerQGEqn::SolveQ( Mat M, Vec f, Vec b, Vector* rhs, Field* qi ) {
	Vec 	x;
	KSP 	ksp;
	Vector*	sol;

	VecDuplicate( f, &x );
	sol = new Vector( "q1-sol", qi, x, NULL, 0 );
	VecZeroEntries( x );
	VecZeroEntries( f );
	rhs->Assemble();
	VecAXPY( f, 1.0, b );
	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, M, M, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "qi_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	KSPDestroy( ksp );
	VecDestroy( x );
	delete sol;
}

/* RHS operator */
DissipationRHS::DissipationRHS( string _name, Mesh* _mesh, double _constant, Field* _field, Field* _tau, Field* _omegaHat, double _rF, double _Ei, Field* _tau_e, double _nu ) : 
RHSOp( _name, _mesh, _constant, _field ) {
	tau		= _tau;
	omegaHat	= _omegaHat;
	rF		= _rF;
	Ei		= _Ei;
	tau_e		= _tau_e;
	nu		= _nu;

	gW = new double*[1];
	gW[0] = new double[2];
}

DissipationRHS::~DissipationRHS() {
	delete[] gW[0];
	delete[] gW;
}

void DissipationRHS::AssembleElement( int el_i, double* rhs ) {
	double	*Ni, **GNx, *coord, detJac, weight, c, q, t, w, te, gCoord[2];

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord	= mesh->el->quadPts[pt_i]->coord;
		weight	= mesh->el->quadPts[pt_i]->weight;
		Ni 	= mesh->el->ShapeFuncs( pt_i );
		GNx	= mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );

		c = detJac*weight*constant;

		mesh->LocalToGlobal( coord, el_i, gCoord );
		field->InterpLocal( el_i, coord, &q );
		tau->InterpLocal( el_i, coord, &t );
		omegaHat->InterpLocal( el_i, coord, &w );
		omegaHat->InterpDerivsGlobal( gCoord, gW );
		tau_e->InterpLocal( el_i, coord, &te );

		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			rhs[node_i] += detJac*weight*q*Ni[node_i];
			rhs[node_i] += c*rF*t*Ni[node_i];
			rhs[node_i] -= c*Ei*w*Ni[node_i];
			rhs[node_i] -= c*rF*te*Ni[node_i];
			rhs[node_i] -= c*nu*( GNx[0][node_i]*gW[0][0] + GNx[1][node_i]*gW[0][1] );
		}
	}
}
