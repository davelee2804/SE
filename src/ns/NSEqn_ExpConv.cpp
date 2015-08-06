#include <iostream>
#include <string>

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
#include "SVV.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "Advector.h"
#include "OIFS.h"
#include "NSEqn_ExpConv.h"

using namespace std;
using std::string;

/* solves the 2D incompressible Navier-Stokes equations
   reference:
	Xu, J. D. Xiu & G. E. Karniadakis (2002) "A Semi-Lagrangian Method for Turbulence 
	Simulations Using Mixed Spectral Discretisations" Journal of Scientific Computing
	vol. 17, 585 - 597 */

NSEqn_ExpConv::NSEqn_ExpConv( Field* velocity, Field* pressure, double _dt ) {
	int 		vSize 	= 2*velocity->mesh->nVertsTotal - velocity->bcs->size[0] - velocity->bcs->size[1];
	int 		pSize 	= pressure->mesh->nVertsTotal - pressure->bcs->size[0];
	Laplacian* 	lapMat;
	Operator**	ops;
	Matrix*		opMat;
	Vector*		solnVec;
	Vector*		rhsVec;

	dt = _dt;
	firstOrderInit = false;

	velStep1 = new Field( "velStep1", velocity->mesh, 2, NULL );
	velStep2 = new Field( "velStep2", velocity->mesh, 2, NULL );

	InitMats( vSize, pSize, velocity->mesh->el->nNodes, pressure->mesh->el->nNodes );

	/* assemble the pressure equation laplacian matrix */
	InitVecs( pSize );
	VecDuplicate( b, &bp );
	VecZeroEntries( bp );
	solnVec = new Vector( "xp", pressure, x, NULL, 0 );
	rhsVec = new Vector( "bp", pressure, bp, NULL, 0 );
	lapMat = new Laplacian( "lapMat", pressure, pressure, -dt );
	ops = new Operator*[1];
	ops[0] = lapMat;
	opMat = new Matrix( "LaplacianOp", Ap, solnVec, solnVec, rhsVec, ops, 1 );
	opMat->Assemble();
	MatAssemblyBegin( Ap, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( Ap, MAT_FINAL_ASSEMBLY );
	delete opMat;
	delete solnVec;
	delete rhsVec;
	VecDestroy( x );
	VecDestroy( b );
}

NSEqn_ExpConv::~NSEqn_ExpConv() {
	delete velStep1;
	delete velStep2;

	MatDestroy( Av );
	MatDestroy( Ap );
	VecDestroy( bv );
	VecDestroy( bp );
}

void NSEqn_ExpConv::Solve( Field* velocity, Field* pressure, Field* prevVel, bool doVelMat, double nu, int svvCutoff ) {
	if( prevVel ) {
		SolveSecondOrder( velocity, pressure, prevVel, doVelMat, nu, svvCutoff );
	}
	else {
		SolveFirstOrder( velocity, pressure, nu, svvCutoff );
	}
}

void NSEqn_ExpConv::SolveFirstOrder( Field* velocity, Field* pressure, double nu, int svvCutoff ) {
	Advector* 	adv;
	int		node_i;
	DivVelRHS*	divVelRHS;
	RHSOp**		rhs;
	Vector*		solnVec;
	Vector*		rhsVec;
	KSP		ksp;
	MatNullSpace	null;
	double		**gP;
	FieldRHS*	velRHS;
	SVV*		svvVisc;
	Laplacian2D*	lapVisc;
	MassMatrix*	massMat;
	Operator**	ops;
	Matrix*		opMat;
	int		nOps	= 1;
	int		op_i	= 0;

	firstOrderInit = true;

	cout << "advective step...\n";
	adv = new Advector( velocity, velocity );
	adv->Advect( dt );
	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		velStep1->vals[node_i][0] = adv->fieldSL->vals[node_i][0];
		velStep1->vals[node_i][1] = adv->fieldSL->vals[node_i][1];
	}
	delete adv;

	cout << "pressure step...\n";
	InitVecs( pressure->mesh->nVertsTotal - pressure->bcs->size[0] );
	solnVec = new Vector( "x-pressureStep", pressure, x, NULL, 0 );

	divVelRHS = new DivVelRHS( "divVelRHS", pressure->mesh, 1.0, velStep1 );
	rhs = new RHSOp*[1];
	rhs[0] = divVelRHS;
	rhsVec = new Vector( "b-pressureStep", pressure, b, rhs, 1 );
	rhsVec->Assemble();
	VecAXPY( b, 1.0, bp );
	VecAssemblyBegin( b );
	VecAssemblyEnd( b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, Ap, Ap, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "p0_" );
	KSPSetFromOptions( ksp );
	MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
	KSPSetNullSpace( ksp, null );
	KSPSolve( ksp, b, x );
	solnVec->UpdateField();
	MatNullSpaceDestroy( null );
	KSPDestroy( ksp );
	VecDestroy( x );
	VecDestroy( b );
	delete solnVec;
	delete rhsVec;

	gP = new double*[1];
	gP[0] = new double[2];
	for( node_i = 0; node_i < velStep2->mesh->nVertsTotal; node_i++ ) {
		pressure->InterpDerivsGlobal( velStep2->mesh->verts[node_i], gP );
		velStep2->vals[node_i][0] = velStep1->vals[node_i][0] - dt*gP[0][0];
		velStep2->vals[node_i][1] = velStep1->vals[node_i][1] - dt*gP[0][1];
	}
	delete[] gP[0];
	delete[] gP;

	cout << "viscous step...\n";
	InitVecs( 2*velocity->mesh->nVertsTotal - velocity->bcs->size[0] - velocity->bcs->size[1] );
	VecDuplicate( b, &bv );
	VecZeroEntries( bv );
	solnVec = new Vector( "x-velocityStep", velocity, x, NULL, 0 );
	/* rhs setup */
	velRHS = new FieldRHS( "velRHS", velocity->mesh, 1.0, velStep2 );
	rhs = new RHSOp*[1];
	rhs[0] = velRHS;
	rhsVec = new Vector( "b-velocityStep", velocity, b, rhs, 1 );
	/* matrix setup */
	massMat = new MassMatrix( "massMat", velocity, velocity, 1.0 );
	if( nu > 0 ) {
		nOps++;
		lapVisc = new Laplacian2D( "lapVisc", velocity, velocity, dt*nu );
	}
	if( svvCutoff > 0 ) {
		nOps++;
		svvVisc = new SVV( "svvVisc", velocity, velocity, dt, velocity->mesh->el->N/2 );
	}
	ops = new Operator*[nOps];
	ops[op_i] = massMat;
	if( nu > 0 ) {
		op_i++;
		ops[op_i] = lapVisc;
	}
	if( svvCutoff > 0 ) {
		op_i++;
		ops[op_i] = svvVisc;
	}
	opMat = new Matrix( "HelmholtzOp", Av, solnVec, solnVec, rhsVec, ops, nOps );
	/* matrix assembly */
	opMat->Assemble();
	MatAssemblyBegin( Av, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( Av, MAT_FINAL_ASSEMBLY );
	/* copy the bcs */
	VecAXPY( bv, 1.0, b );
	rhsVec->Assemble();
	VecAssemblyBegin( b );
	VecAssemblyEnd( b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, Av, Av, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "v0_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, b, x );
	solnVec->UpdateField();
	KSPDestroy( ksp );
	VecDestroy( x );
	VecDestroy( b );
	delete opMat;
	delete solnVec;
	delete rhsVec;
}

void NSEqn_ExpConv::SolveSecondOrder( Field* velocity, Field* pressure, Field* prevVel, bool doVelMat, double nu, int svvCutoff ) {
	ConvectionMatrix*	convMat;
	OIFS*			expConv;
	int			node_i;
	DivVelRHS*		divVelRHS;
	RHSOp**			rhs;
	Vector*			solnVec;
	Vector*			rhsVec;
	KSP			ksp;
	MatNullSpace		null;
	double			**gP;
	FieldRHS*		velRHS;
	SVV*			svvVisc;
	Laplacian2D*		lapVisc;
	MassMatrix*		massMat;
	Operator**		ops;
	Matrix*			opMat	= NULL;
	int			nOps	= 1;
	int			op_i	= 0;

	cout << "advective step...\n";

	convMat = new ConvectionMatrix( "convMat", velocity, velocity, velocity );
	expConv = new OIFS( convMat, velocity, prevVel, velocity, prevVel, dt, 1 );
	expConv->Solve();
	velStep1->Copy( expConv->phiTilde );
	delete expConv;

	cout << "pressure step...\n";
	InitVecs( pressure->mesh->nVertsTotal - pressure->bcs->size[0] );
	solnVec = new Vector( "x-pressureStep", pressure, x, NULL, 0 );

	divVelRHS = new DivVelRHS( "divVelRHS", pressure->mesh, 1.0, velStep1 );
	rhs = new RHSOp*[1];
	rhs[0] = divVelRHS;
	rhsVec = new Vector( "b-pressureStep", pressure, b, rhs, 1 );
	rhsVec->Assemble();
	VecAXPY( b, 1.0, bp );
	VecAssemblyBegin( b );
	VecAssemblyEnd( b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, Ap, Ap, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "p1_" );
	KSPSetFromOptions( ksp );
	MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
	KSPSetNullSpace( ksp, null );
	KSPSolve( ksp, b, x );
	solnVec->UpdateField();
	MatNullSpaceDestroy( null );
	KSPDestroy( ksp );
	VecDestroy( x );
	VecDestroy( b );
	delete solnVec;
	delete rhsVec;

	gP = new double*[1];
	gP[0] = new double[2];
	for( node_i = 0; node_i < velStep2->mesh->nVertsTotal; node_i++ ) {
		pressure->InterpDerivsGlobal( velStep2->mesh->verts[node_i], gP );
		velStep2->vals[node_i][0] = velStep1->vals[node_i][0] - dt*gP[0][0];
		velStep2->vals[node_i][1] = velStep1->vals[node_i][1] - dt*gP[0][1];
	}
	delete[] gP[0];
	delete[] gP;

	cout << "viscous step...\n";
	InitVecs( 2*velocity->mesh->nVertsTotal - velocity->bcs->size[0] - velocity->bcs->size[1] );
	if( doVelMat ) {
		if( firstOrderInit ) {
			VecDestroy( bv );
		}
		VecDuplicate( b, &bv );
		VecZeroEntries( bv );
		MatZeroEntries( Av );
	}
	solnVec = new Vector( "x-velocityStep", velocity, x, NULL, 0 );
	velRHS = new FieldRHS( "velRHS", velocity->mesh, 1.0, velStep2 );
	rhs = new RHSOp*[1];
	rhs[0] = velRHS;
	rhsVec = new Vector( "b-velocityStep", velocity, b, rhs, 1 );
	if( doVelMat ) {
		massMat = new MassMatrix( "massMat", velocity, velocity, 1.5 );
		if( nu > 0 ) {
			nOps++;
			lapVisc = new Laplacian2D( "lapVisc", velocity, velocity, dt*nu );
		}
		if( svvCutoff > 0 ) {
			nOps++;
			svvVisc = new SVV( "svvVisc", velocity, velocity, dt, velocity->mesh->el->N/2 );
		}
		ops = new Operator*[nOps];
		ops[op_i] = massMat;
		if( nu > 0 ) {
			op_i++;
			ops[op_i] = lapVisc;
		}
		if( svvCutoff > 0 ) {
			op_i++;
			ops[op_i] = svvVisc;
		}
		opMat = new Matrix( "HelmholtzOp", Av, solnVec, solnVec, rhsVec, ops, nOps );
		opMat->Assemble();
		MatAssemblyBegin( Av, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd( Av, MAT_FINAL_ASSEMBLY );
		VecAXPY( bv, 1.0, b );
	}
	rhsVec->Assemble();
	/* copy the bcs */
	if( !doVelMat ) {
		VecAXPY( b, 1.0, bv );
	}
	VecAssemblyBegin( b );
	VecAssemblyEnd( b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, Av, Av, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "v1_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, b, x );
	solnVec->UpdateField();
	KSPDestroy( ksp );
	VecDestroy( x );
	VecDestroy( b );
	if( doVelMat ) {
		delete opMat;
	}
	delete solnVec;
	delete rhsVec;
}

void NSEqn_ExpConv::InitMats( int vSize, int pSize, int nVelNodes, int nPresNodes ) {
	MatCreate( MPI_COMM_WORLD, &Av );
	MatSetSizes( Av, vSize, vSize, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( Av, MATSEQAIJ );
	MatSeqAIJSetPreallocation( Av, 4*2*nVelNodes, PETSC_NULL );
	MatZeroEntries( Av );

	MatCreate( MPI_COMM_WORLD, &Ap );
	MatSetSizes( Ap, pSize, pSize, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( Ap, MATSEQAIJ );
	MatSeqAIJSetPreallocation( Ap, 4*nPresNodes, PETSC_NULL );
	MatZeroEntries( Ap );
}

void NSEqn_ExpConv::InitVecs( int size ) {
	VecCreate( MPI_COMM_WORLD, &b );
	VecSetSizes( b, size, PETSC_DETERMINE );
	VecSetType( b, VECSEQ );
	VecSetOption( b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
	VecZeroEntries( b );

	VecDuplicate( b, &x );
	VecZeroEntries( x );
}
