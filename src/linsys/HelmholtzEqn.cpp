#include <cstdlib>
#include <string>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "HelmholtzEqn.h"

using namespace std;
using std::string;

/* solves the Helmholtz equation for the unknown vector x, given a forcing term as represented by a scalar 
 * field. assumes that the vector x has already been build and had its boundary conditions applied, and that
 * the field has been initialised. */

HelmholtzEqn::HelmholtzEqn( Field* _solField, Field* _rhsField, double _c, double _nu, double _dt, int _svvCutoff ) {
	MassMatrix*	masOp;
	Laplacian*	lapOp;
	SVV*		svvOp	= NULL;
	FieldRHS*	rhsOp	= NULL;
	int 		size	= 0;
	Operator**	matOps;
	RHSOp**		rhsOps	= NULL;
	int		nOps	= (_svvCutoff) ? 3 : 2;
	int		nRHSOps = (_rhsField) ? 1 : 0;

	matOps = new Operator*[nOps];
	if( nRHSOps ) { rhsOps = new RHSOp*[1]; }

	solField = _solField;
	rhsField = _rhsField;
	c  = _c;
	nu = _nu;
	dt = _dt;
	svvCutoff = _svvCutoff;

	assembled = false;

	for( int dof_i = 0; dof_i < solField->nDofs; dof_i++ ) {
		size += solField->mesh->nVertsTotal - solField->bcs->size[dof_i];
	}

	/* setup the petsc objects */
	MatCreate( MPI_COMM_WORLD, &A );
	MatSetSizes( A, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( A, MATSEQAIJ );
	MatSetFromOptions( A );
	MatSeqAIJSetPreallocation( A, 4*solField->nDofs*solField->nDofs*solField->mesh->el->nNodes, PETSC_NULL );
	MatZeroEntries( A );
	
	VecCreate( MPI_COMM_WORLD, &x );
	VecSetSizes( x, size, PETSC_DETERMINE );
	VecSetType( x, VECSEQ );
	VecSetFromOptions( x );
	VecSetOption( x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
	VecZeroEntries( x );

	VecDuplicate( x, &f );
	VecZeroEntries( f );

	VecDuplicate( x, &b );
	VecZeroEntries( b );

	lapOp = new Laplacian( "lapOp", solField, solField, nu*dt );
	masOp = new MassMatrix( "massOp", solField, solField, c );
	if( svvCutoff ) { svvOp = new SVV( "svvOp", solField, solField, dt, svvCutoff ); }
	if( nRHSOps )   { rhsOp = new FieldRHS( "rhs", solField->mesh, 1.0, rhsField );  }

	matOps[0] = masOp;
	matOps[1] = lapOp;
	if( svvCutoff ) { matOps[2] = svvOp; }
	if( nRHSOps )   { rhsOps[0] = rhsOp; }

	sol = new Vector( "x", solField, x, NULL, 0 );
	rhs = new Vector( "f", solField, f, rhsOps, nRHSOps );
	op  = new Matrix( "A", A, sol, sol, rhs, matOps, nOps );
}

HelmholtzEqn::~HelmholtzEqn() {
	MatDestroy( &A );
	VecDestroy( &f );
	VecDestroy( &x );
	VecDestroy( &b );
}

void HelmholtzEqn::Assemble() {
	MatZeroEntries( A );
	VecZeroEntries( f );
	op->Assemble();
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );

	/* copy the bcs */
	VecZeroEntries( b );
	VecCopy( f, b );

	assembled = true;
}

void HelmholtzEqn::Solve( std::string name, bool removeNullSpace ) {
	KSP 		ksp;
	MatNullSpace	null;

	if( !assembled ) { Assemble(); }

	VecZeroEntries( f );
	rhs->Assemble();
	VecAssemblyBegin( f );
	VecAssemblyEnd( f );

	/* copy the bcs */
	VecAXPY( f, 1.0, b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	//KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
	KSPSetOperators( ksp, A, A );
	KSPSetOptionsPrefix( ksp, name.c_str() );
	KSPSetFromOptions( ksp );
	if( removeNullSpace ) {
		MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
		KSPSetNullSpace( ksp, null );
	}
	KSPSolve( ksp, f, x );
	KSPDestroy( &ksp );
	if( removeNullSpace ) {
		MatNullSpaceDestroy( &null );
	}

	/* copy the solution vector values onto the corresponding field */
	sol->UpdateField();
}
