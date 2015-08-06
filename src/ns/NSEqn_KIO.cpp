#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "NSEqn_KIO.h"

using namespace std;
using std::string;

NSEqn_KIO::NSEqn_KIO( Field* _velocity, Field* _pressure, double _dt, double _nu ) {
	velocity = _velocity;
	pressure = _pressure;
	dt       = _dt;
	nu       = _nu;

	M_u 	= NULL;
	M_p 	= NULL;
	b_u 	= NULL;
}

NSEqn_KIO::~NSEqn_KIO() {
	MatDestroy( &M_u );
	MatDestroy( &M_p );
	VecDestroy( &b_u );
}

void NSEqn_KIO::Setup( int order ) {
	MassMatrix*		massMat;
	Laplacian*		viscMat;
	Laplacian*		lapMat;
	Operator**		ops;
	Vector*			uSol;
	Vector*			pSol;
	Vector*			uRhs;
	Matrix*			uMat;
	Matrix*			pMat;
	Vec			u;
	Vec			p;
	double			alpha		= ( order == 2 ) ? 1.5 : 1.0;

	if( M_u ) { MatDestroy( &M_u ); }
	if( M_p ) { MatDestroy( &M_p ); }
	if( b_u ) { VecDestroy( &b_u ); }

	CreateMatrix( velocity, velocity, &M_u );
	CreateMatrix( pressure, pressure, &M_p );

	CreateVector( velocity, &u );
	CreateVector( pressure, &p );
	CreateVector( velocity, &b_u );

	uSol = new Vector( "u-sol", velocity, u, NULL, 0 );
	pSol = new Vector( "p-sol", pressure, p, NULL, 0 );
	uRhs = new Vector( "u-rhs", velocity, b_u, NULL, 0 );

	/* assemble the velocity step matrix */
	massMat = new MassMatrix( "u-mass-mat", velocity, velocity, alpha );
	viscMat = new Laplacian( "u-visc-mat", velocity, velocity, dt*nu );
	ops = new Operator*[2];
	ops[0] = massMat;
	ops[1] = viscMat;
	uMat = new Matrix( "M_u", M_u, uSol, uSol, uRhs, ops, 2 );
	uMat->Assemble();
	MatAssemblyBegin( M_u, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( M_u, MAT_FINAL_ASSEMBLY );

	/* assemble the pressure step matrix */
	lapMat = new Laplacian( "p-lap-mat", pressure, pressure, -dt/alpha );
	ops = new Operator*[1];
	ops[0] = lapMat;
	pMat = new Matrix( "M_p", M_p, pSol, pSol, NULL, ops, 1 );
	pMat->Assemble();
	MatAssemblyBegin( M_p, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( M_p, MAT_FINAL_ASSEMBLY );

	delete uSol;
	delete pSol;
	delete uRhs;
	delete uMat;
	delete pMat;
	VecDestroy( &u );
	VecDestroy( &p );
}

void NSEqn_KIO::Solve( Field* velPrev ) {
	Field*			velHat		= new Field( "vel-hat", velocity->mesh, 2, velocity->bcs );
	Field*			presDelta	= new Field( "pres-delta", pressure->mesh, 1, pressure->bcs );
	FieldRHS*		velCurrRHS;
	FieldRHS*		velPrevRHS	= NULL;
	ConvectionRHS*		convCurrRHS;
	ConvectionRHS*		convPrevRHS	= NULL;
	GradFieldRHS*		gradPresRHS;
	DivVelRHS*		divVelRHS;
	RHSOp**			ops;
	Vector*			uSol;
	Vector*			pSol;
	Vector*			uRhs;
	Vector*			pRhs;
	Vec			u;
	Vec			p;
	Vec			fu;
	Vec			fp;
	double			alpha		= ( velPrev ) ? 1.5 : 1.0;
	double			delta		= ( velPrev ) ? 2.0 : 1.0;
	int			nVelOps		= ( velPrev ) ? 5 : 3;
	KSP			ksp;
	double			**dp;
	MatNullSpace		null;

	CreateVector( velocity, &u );
	CreateVector( pressure, &p );
	CreateVector( velocity, &fu );
	CreateVector( pressure, &fp );

	uSol = new Vector( "u-sol", velHat, u, NULL, 0 );
	pSol = new Vector( "p-sol", presDelta, p, NULL, 0 );

	/* assemble the velocity equation rhs vector */
	ops = new RHSOp*[nVelOps];
	velCurrRHS  = new FieldRHS( "vel-curr-rhs", velocity->mesh, delta, velocity );
	convCurrRHS = new ConvectionRHS( "conv-curr-rhs", velocity->mesh, -delta*dt, velocity );
	gradPresRHS = new GradFieldRHS( "grad-pres-rhs", velocity->mesh, -dt, pressure );
	ops[0] = velCurrRHS;
	ops[1] = convCurrRHS;
	ops[2] = gradPresRHS;
	if( velPrev ) { 
		velPrevRHS  = new FieldRHS( "vel-prev-rhs", velocity->mesh, -0.5, velPrev ); 
		convPrevRHS = new ConvectionRHS( "conv-prev-rhs", velocity->mesh, dt, velPrev );
		ops[3] = velPrevRHS;
		ops[4] = convPrevRHS;
	}
	uRhs = new Vector( "v-rhs", velHat, fu, ops, nVelOps );
	uRhs->Assemble();

	VecAXPY( fu, 1.0, b_u ); /* ??? for periodic bcs? */

	/* solve the velocity equation */
	KSPCreate( MPI_COMM_WORLD, &ksp );
	//KSPSetOperators( ksp, M_u, M_u, SAME_NONZERO_PATTERN );
	KSPSetOperators( ksp, M_u, M_u );
	KSPSetOptionsPrefix( ksp, "helm_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, fu, u );
	KSPDestroy( &ksp );
	uSol->UpdateField();

	/* assemble the pressure equation rhs vector */
	divVelRHS = new DivVelRHS( "div-height-vel-rhs", pressure->mesh, 1.0, velHat );
	ops = new RHSOp*[1];
	ops[0] = divVelRHS;
	pRhs = new Vector( "p-rhs", presDelta, fp, ops, 1 );
	pRhs->Assemble();

	/* solve the pressure equation */
	KSPCreate( MPI_COMM_WORLD, &ksp );
	//KSPSetOperators( ksp, M_p, M_p, SAME_NONZERO_PATTERN );
	KSPSetOperators( ksp, M_p, M_p );
	KSPSetOptionsPrefix( ksp, "poisson_" );
	KSPSetFromOptions( ksp );
	MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
	KSPSetNullSpace( ksp, null );
	KSPSolve( ksp, fp, p );
	KSPDestroy( &ksp );
	MatNullSpaceDestroy( &null );
	pSol->UpdateField();

	/* update the velocity */
	dp = new double*[1];
	dp[0] = new double[2];
	for( int node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		presDelta->InterpDerivsGlobal( velocity->mesh->verts[node_i], dp );
		velocity->vals[node_i][0] = velHat->vals[node_i][0] - (dt/alpha)*dp[0][0];
		velocity->vals[node_i][1] = velHat->vals[node_i][1] - (dt/alpha)*dp[0][1];
	}
	delete[] dp[0];
	delete[] dp;

	for( int node_i = 0; node_i < pressure->mesh->nVertsTotal; node_i++ ) {
		pressure->vals[node_i][0] += presDelta->vals[node_i][0];
	}

	delete velHat;
	delete presDelta;

	delete uSol;
	delete pSol;
	delete uRhs;
	delete pRhs;

	VecDestroy( &u );
	VecDestroy( &p );
	VecDestroy( &fu );
	VecDestroy( &fp );
}
