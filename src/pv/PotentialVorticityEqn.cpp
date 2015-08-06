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
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "Advector.h"
#include "PotentialVorticityEqn.h"

using namespace std;
using std::string;

PotentialVorticityEqn::PotentialVorticityEqn( Field* _omega, Field* _psi, Field* _velocity, Field* _prevOmega, Field* _prevVel ) {
	MassMatrix* omegaMassMat;
	Laplacian* laplacian;
	MassMatrix* psiMassMat;
	FieldRHS* omega1RHS;
	FieldRHS* omega2RHS = NULL;
	Operator** ops;
	RHSOp** rhs;
	static double rhsFac;
	static int size;

	omega = _omega;
	psi = _psi;
	velocity = _velocity;
	prevOmega = _prevOmega;
	prevVel = _prevVel;

	size = omega->mesh->nVertsTotal - omega->bcs->size[0] + psi->mesh->nVertsTotal - psi->bcs->size[0];

	/* set up the petsc objects */
	MatCreate( MPI_COMM_WORLD, &A );
	MatSetSizes( A, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( A, MATSEQAIJ );
	MatSetFromOptions( A );
	MatSeqAIJSetPreallocation( A, 4*4*omega->mesh->el->nNodes, PETSC_NULL );
	MatZeroEntries( A );

	VecCreate( MPI_COMM_WORLD, &x );
	VecSetSizes( x, size, PETSC_DETERMINE );
	VecSetType( x, VECSEQ );
	VecSetFromOptions( x );
	VecSetOption( x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
	VecZeroEntries( x );

	VecCreate( MPI_COMM_WORLD, &b );
	VecSetSizes( b, size, PETSC_DETERMINE );
	VecSetType( b, VECSEQ );
	VecSetFromOptions( b );
	VecSetOption( b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
	VecZeroEntries( b );

	if( prevOmega && prevVel ) {
		advector = new Advector( omega, velocity, prevOmega, prevVel );
		rhsFac = 4.0/3.0;
	}
	else {
		advector = new Advector( omega, velocity );
		rhsFac = 1.0;
	}

	omega1RHS = new FieldRHS( "omega1RHS", omega->mesh, rhsFac, advector->fieldSL );
	if( prevOmega && prevVel ) {
		omega2RHS = new FieldRHS( "omega2RHS", omega->mesh, -1.0/3.0, advector->fieldSLMinusOne );
	}

	x0 = new Vector( "xOmega", omega, x, NULL, 0 );
	x1 = new Vector( "xPsi", psi, x, NULL, 0 );
	
	if( prevOmega && prevVel ) {
		rhs = new RHSOp*[2];
		rhs[0] = omega1RHS;
		rhs[1] = omega2RHS;
		b0 = new Vector( "fOmega", omega, b, rhs, 2 );
	}
	else {
		rhs = new RHSOp*[1];
		rhs[0] = omega1RHS;
		b0 = new Vector( "fOmega", omega, b, rhs, 1 );
	}
	b1 = new Vector( "fPsi", psi, b, NULL, 0 );

	omegaMassMat = new MassMatrix( "omegaMassMat", omega, omega, 1.0 );
	psiMassMat = new MassMatrix( "psiMassMat", psi, omega, 1.0 );
	laplacian = new Laplacian( "laplacian", psi, psi, 1.0 );

	ops = new Operator*[1];
	ops[0] = omegaMassMat;
	A00 = new Matrix( "A00", A, x0, x0, b0, ops, 1 );

	ops = new Operator*[1];
	ops[0] = psiMassMat;
	A10 = new Matrix( "A10", A, x1, x0, b1, ops, 1 );

	ops = new Operator*[1];
	ops[0] = laplacian;
	A11 = new Matrix( "A11", A, x1, x1, b1, ops, 1 );
}

PotentialVorticityEqn::~PotentialVorticityEqn() {
	VecDestroy( x );
	VecDestroy( b );
	MatDestroy( A );

	delete advector;

	delete x0;
	delete x1;

	delete b0;
	delete b1;

	delete A00;
	delete A10;
	delete A11;
}

void PotentialVorticityEqn::Solve( double dt ) {
	KSP ksp;
	PC pc;
	IS omegaIS, psiIS;

	advector->Advect( dt );

	A00->Assemble();
	A10->Assemble();
	A11->Assemble();
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );

	b0->Assemble();
	VecAssemblyBegin( b );
	VecAssemblyEnd( b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
	KSPGetPC( ksp, &pc );
	PCSetType( pc, PCFIELDSPLIT );
	PCFieldSplitSetType( pc, PC_COMPOSITE_MULTIPLICATIVE );
	ISCreateStride( MPI_COMM_WORLD, x0->vecSize, 0, 1, &omegaIS );
	ISCreateStride( MPI_COMM_WORLD, x1->vecSize, x0->vecSize, 1, &psiIS );
	PCFieldSplitSetIS( pc, omegaIS );
	PCFieldSplitSetIS( pc, psiIS );

	KSPSetFromOptions( ksp );
	KSPSolve( ksp, b, x );
	KSPDestroy( ksp );

	ISDestroy( omegaIS );
	ISDestroy( psiIS );

	x0->UpdateField();
	x1->UpdateField();
}
