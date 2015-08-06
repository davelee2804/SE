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
#include "StokesEqn.h"

using namespace::std;
using std::string;

StokesEqn::StokesEqn( Field* _velocity, Field* _pressure, Field* _forcing ) {
	int vSize, pSize, size;
	Laplacian* 	laplacian;
	Gradient*	gradient;
	Divergence*	divergence;
	StokesPC*	stokesPC;
	FieldRHS*	fieldRHS;
	Operator**	matOps;
	RHSOp**		rhsOps;

	velocity = _velocity;
	pressure = _pressure;
	forcing = _forcing;

	vSize = 2*velocity->mesh->nVertsTotal - velocity->bcs->size[0] - velocity->bcs->size[1];
	pSize = pressure->mesh->nVertsTotal;
	size = vSize + pSize;

	/* set up the petsc objects */
	MatCreate( MPI_COMM_WORLD, &A );
	MatSetSizes( A, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( A, MATSEQAIJ );
	MatSetFromOptions( A );
	MatSeqAIJSetPreallocation( A, 4*4*2*velocity->mesh->el->nNodes, PETSC_NULL );
	MatZeroEntries( A );

	MatCreate( MPI_COMM_WORLD, &P );
	MatSetSizes( P, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( P, MATSEQAIJ );
	MatSetFromOptions( P );
	MatSeqAIJSetPreallocation( P, 4*4*2*velocity->mesh->el->nNodes, PETSC_NULL );
	MatZeroEntries( P );

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

	laplacian = new Laplacian( "laplacian", velocity, velocity, 1.0 );
	gradient = new Gradient( "gradient", velocity, pressure, 1.0 );
	divergence = new Divergence( "divergence", pressure, velocity, 1.0 );
	stokesPC = new StokesPC( "stokesPC", pressure, pressure, 1.0, velocity, pressure );

	fieldRHS = new FieldRHS( "stokesForcing", velocity->mesh, 1.0, forcing );

	v = new Vector( "v", velocity, x, NULL, 0 );
	p = new Vector( "p", pressure, x, NULL, 0 );

	rhsOps = new RHSOp*[1];
	rhsOps[0] = fieldRHS;
	f = new Vector( "f", velocity, b, rhsOps, 1 );
	h = new Vector( "h", pressure, b, NULL, 0 );

	matOps = new Operator*[1];
	matOps[0] = laplacian;
	K = new Matrix( "K", A, v, v, f, matOps, 1 );

	matOps = new Operator*[1];
	matOps[0] = gradient;
	G = new Matrix( "G", A, v, p, f, matOps, 1 );

	matOps = new Operator*[1];
	matOps[0] = divergence;
	D = new Matrix( "D", A, p, v, h, matOps, 1 );

	matOps = new Operator*[1];
	matOps[0] = stokesPC;
	C = new Matrix( "C", P, p, p, NULL, matOps, 1 );
}

StokesEqn::~StokesEqn() {
	MatDestroy( &A );
	MatDestroy( &P );
	VecDestroy( &x );
	VecDestroy( &b );

	delete K;
	delete G;
	delete D;
	delete C;

	delete v;
	delete p;

	delete f;
	delete h;
}

void StokesEqn::Solve() {
	KSP ksp;
	PC pc;
	MatNullSpace null;
	IS vis, pis;
	KSP* subksp;
	PetscInt n;

	/* assemble the stokes problem */
	K->Assemble();
	G->Assemble();
	D->Assemble();
	/* assemble the pressure-pressure preconditioner block */
	C->Assemble();
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
	/* and copy across into the preconditioner matrix */
	MatAXPY( P, 1.0, A, DIFFERENT_NONZERO_PATTERN );
	MatAssemblyBegin( P, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( P, MAT_FINAL_ASSEMBLY );
	/* assemble the rhs forcing */
	f->Assemble();
	VecAssemblyBegin( b );
	VecAssemblyEnd( b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	//KSPSetOperators( ksp, A, P, DIFFERENT_NONZERO_PATTERN );
	KSPSetOperators( ksp, A, P );
	KSPGetPC( ksp, &pc );
	PCSetType( pc, PCFIELDSPLIT );
	PCFieldSplitSetType( pc, PC_COMPOSITE_MULTIPLICATIVE );
	ISCreateStride( MPI_COMM_WORLD, v->vecSize, 0, 1, &vis );
	ISCreateStride( MPI_COMM_WORLD, p->vecSize, v->vecSize, 1, &pis );
	//PCFieldSplitSetIS( pc, vis );
	//PCFieldSplitSetIS( pc, pis );
	PCFieldSplitSetIS( pc, "0", vis );
	PCFieldSplitSetIS( pc, "1", pis );

	/* remove the pressure null space */
	PCFieldSplitGetSubKSP( pc, &n, &subksp );
	MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
	KSPSetNullSpace( subksp[1], null );

	KSPSetFromOptions( ksp );
	KSPSolve( ksp, b, x );
	KSPDestroy( &ksp );

	MatNullSpaceDestroy( &null );
	ISDestroy( &vis );
	ISDestroy( &pis );

	v->UpdateField();
	p->UpdateField();
}
