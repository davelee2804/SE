#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SchurComplement.h"

/*
Solves a mixed (ie: velocity-pressure) linear system using the Uzawa method:
	[K G][u] = [f]
	[D C][p] = [h]

	=> [DK^{-1}G - C][p] = [DK^{-1}][f] - [h]
	=>               [u] = [K^{-1}][f] - [K^{-1}G][p]

requires an operator in order to assemble the matrix [DK^{-1}G - C]
*/

using namespace std;
using std::string;

SchurComplement::SchurComplement( Mat _Kinv, Mat _G, Mat _D, Mat _C ) {
	Kinv = _Kinv;
	G    = _G;
	D    = _D;
	C    = _C;

	MatMatMult( D, Kinv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DKinv );
	MatMatMult( DKinv, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &DKinvG );
	if( C ) { MatAXPY( DKinvG, -1.0, C, SUBSET_NONZERO_PATTERN ); } /* ??? */
}

SchurComplement::~SchurComplement() {
	MatDestroy( &DKinv );
	MatDestroy( &DKinvG );
}

void SchurComplement::Solve( string prefix, Vec v, Vec p, Vec f, Vec h ) {
	Vec DKinvf;
	Vec Kinvf;
	Vec Gp;
	Vec KinvGp;
	KSP ksp;

	VecDuplicate( p, &DKinvf );
	VecDuplicate( v, &Kinvf );
	VecDuplicate( v, &Gp );
	VecDuplicate( v, &KinvGp );

	/* setup rhs */
	MatMult( DKinv, f, DKinvf );
	if( h ) VecAXPY( DKinvf, -1.0, h );

	/* solve for p */
	KSPCreate( MPI_COMM_WORLD, &ksp );
	//KSPSetOperators( ksp, DKinvG, DKinvG, SAME_NONZERO_PATTERN );
	KSPSetOperators( ksp, DKinvG, DKinvG );
	KSPSetOptionsPrefix( ksp, prefix.c_str() );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, DKinvf, p );
	KSPDestroy( &ksp );

	/* update v */
	MatMult( Kinv, f, Kinvf );
	MatMult( G, p, Gp );
	MatMult( Kinv, Gp, KinvGp );
	VecZeroEntries( v );
	VecAXPY( v, +1.0, Kinvf );
	VecAXPY( v, -1.0, KinvGp );

	VecDestroy( &DKinvf );
	VecDestroy( &Kinvf );
	VecDestroy( &Gp );
	VecDestroy( &KinvGp );
}
