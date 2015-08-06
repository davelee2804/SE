#include <cstdlib>
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
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "PoissonEqn.h"

using namespace std;
using std::string;

/* solves poissons eqn for the unknown vector x, given a forcing term as represented by a scalar field.
 * assumes that the vector x has already been build and had its boundary conditions applied, and that
 * the field has been initialised.
 */

PoissonEqn::PoissonEqn( Field* _solField, Field* _rhsField, double _c ) {
	Mesh*		mesh;
	Laplacian* 	laplacian;
	FieldRHS*	poissonRHS;
	int 		size;
	Operator**	ops;
	RHSOp**		rhs;

	ops = new Operator*[1];
	rhs = new RHSOp*[1];

	solField = _solField;
	rhsField = _rhsField;
	c = _c;

	mesh = solField->mesh;
	size = mesh->nVertsTotal - solField->bcs->size[0]; /* scalar field */

	/* setup the petsc objects */
	MatCreate( MPI_COMM_WORLD, &A );
	MatSetSizes( A, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( A, MATSEQAIJ );
	MatSetFromOptions( A );
	MatSeqAIJSetPreallocation( A, 4*solField->mesh->el->nNodes, PETSC_NULL );
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

	laplacian = new Laplacian( "laplacian", solField, solField, c );
	poissonRHS = new FieldRHS( "rhs", mesh, 1.0, rhsField );

	ops[0] = laplacian;
	rhs[0] = poissonRHS;

	x0 = new Vector( "x", solField, x, NULL, 0 );
	f0 = new Vector( "f", solField, f, rhs, 1 );
	A00 = new Matrix( "A", A, x0, x0, f0, ops, 1 );

	assembled = false;
}

PoissonEqn::~PoissonEqn() {
	MatDestroy( &A );
	VecDestroy( &f );
	VecDestroy( &x );
}

void PoissonEqn::Assemble() {
	MatZeroEntries( A );
	VecZeroEntries( f );
	A00->Assemble();
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );

	/* copy the bcs */
	VecZeroEntries( b );
	VecCopy( f, b );

	assembled = true;
}

void PoissonEqn::Solve( bool removeNullSpace ) {
	KSP 		ksp;
	MatNullSpace 	null;

	if( !assembled ) {
		Assemble();
	}

	VecZeroEntries( f );
	f0->Assemble();
	VecAssemblyBegin( f );
	VecAssemblyEnd( f );

	VecAXPY( f, 1.0, b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	//KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
	KSPSetOperators( ksp, A, A );
	KSPSetOptionsPrefix( ksp, "poisson_" );
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
	x0->UpdateField();
}
