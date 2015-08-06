#include <cstdlib>
#include <string>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "OIFS.h"

/* Operator-Integration-Factor explicit evaluation of non-linear terms. To be used as the RHS
 * of a linear equation when splitting a system of equations into the linear (implicit) and
 * non-linear (explicit) terms.
 *
 * References:
 *	Maday, Y., A. T. Patera and E. M. Ronquist (1990) "An Operator-Integration-Factor
 *	Splitting Method for Time-Dependent Problems: Applications to Incompressible Fluid
 *	Flow" Journal of Scientific Computing, 5, 263-291
 *
 *	Iskandarani, M., D. B. Haidvogel and J. P. Boyd (1995) "A Staggered Spectral Element
 *	Model with Application to the Oceanic Shallow Water Equations" International Journal
 *	for Numerical Methods in Fluids, 20, 393-414
 *
 *	Timmermans, L. J. P., P. D. Minev and F. N. Van Der Vosse (1996) "An Approximate 
 *	Projection Scheme for Incompressible Flow using Spectral Elements" International
 *	Journal for Numerical Methods in Fluids, 22, 673-688
 */

using namespace std;
using std::string;

OIFS::OIFS( Operator* _op, Field* _phi0, Field* _phi1, Field* _psi0, Field* _psi1, double _dt, int _nCycles ) {
	int 		size, alloc, row, col, nCols;
	MassMatrix* 	massMat;
	Operator** 	ops;
	const 		PetscInt* cols;
	const 		PetscScalar* vals;
	int 		dof_i;

	op = _op;
	phi0 = _phi0; /* the unknown field at the current time step */
	phi1 = _phi1; /* the unknown field at the previous time step */
	psi0 = _psi0; /* the non-linear field at the current time step */
	psi1 = _psi1; /* the non-linear field at the previous time step */
	dt = _dt;
	nCycles = _nCycles;

	/* assumption: the row field and the column field for the non-linear operator are the same, 
	   however the this may differ from the field used to assemble this operator */
	phiTilde = new Field( "phi~", phi0->mesh, phi0->nDofs, phi0->bcs );

	size = phi0->nDofs*phi0->mesh->nVertsTotal;
	for( dof_i = 0; dof_i < phi0->nDofs; dof_i++ ) {
		size -= phi0->bcs->size[dof_i];
	}
	alloc = 4*phi0->nDofs*phi0->mesh->el->nNodes;

	InitMat( &M, size, size, alloc );
	InitMat( &Minv, size, size, alloc );
	InitMat( &A, size, size, alloc );
	InitVec( &x, size );

	vector  = new Vector( "x", phi0, x, NULL, 0 );
	massMat = new MassMatrix( "massMat", phi0, phi0, 1.0 );
	ops     = new Operator*[1];
	ops[0]  = massMat;
	matrix  = new Matrix( "MassMatrix", M, vector, vector, NULL, ops, 1 );
	matrix->Assemble();
	MatAssemblyBegin( M, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( M, MAT_FINAL_ASSEMBLY );

	/* mass matrix inverse assumes shape functions are orthogonal (ie: Legendre), and so M is diagonal */
	for( row = 0; row < size; row++ ) {
		MatGetRow( M, row, &nCols, &cols, &vals );
		for( col = 0; col < nCols; col++ ) {
			if( cols[col] == row ) {
				break;
			}
		}
		MatSetValue( Minv, row, row, 1.0/vals[col], INSERT_VALUES );
		MatRestoreRow( M, row, &nCols, &cols, &vals );
	}
	MatAssemblyBegin( Minv, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( Minv, MAT_FINAL_ASSEMBLY );

	delete matrix;

	ops = new Operator*[1];
	ops[0] = op;
	matrix = new Matrix( "nonLinearOp", A, vector, vector, NULL, ops, 1 );
}

OIFS::~OIFS() {
	MatDestroy( M );
	MatDestroy( Minv );
	MatDestroy( A );
	VecDestroy( x );
	delete vector;
	delete matrix;
	delete phiTilde;
}

void OIFS::Solve() {
	if( phi1 && psi1 ) {
		SolveSecondOrder();
	}
	else {
		SolveFirstOrder();
	}
}

void OIFS::SolveFirstOrder() {
	int cycle_i;
	double ds = dt/nCycles;
	Vec um0, um1, g0, g1, g2, g3, tmp1, tmp2;
	Field* psi = new Field( "psi", psi0->mesh, psi0->nDofs, psi0->bcs );

	cout << "solving for non-linear term " << op->name << ", first order\n";

	op->field = psi;

	VecDuplicate( x, &um0 );
	VecDuplicate( x, &um1 );
	VecDuplicate( x, &g0 );
	VecDuplicate( x, &g1 );
	VecDuplicate( x, &g2 );
	VecDuplicate( x, &g3 );
	VecDuplicate( x, &tmp1 );
	VecDuplicate( x, &tmp2 );

	VecZeroEntries( um0 );
	VecZeroEntries( um1 );
	VecZeroEntries( g0 );
	VecZeroEntries( g1 );
	VecZeroEntries( g2 );
	VecZeroEntries( g3 );

	FieldToVec( phi0, um0 );
	psi->Copy( psi0 );
	MatZeroEntries( A );
	matrix->Assemble();
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
	
	for( cycle_i = 0; cycle_i < nCycles; cycle_i++ ) {
		cout << "\tcycle: " << cycle_i << endl;

		MatMult( A, um0, tmp1 );
		MatMult( Minv, tmp1, g0 );

		VecZeroEntries( tmp1 );
		VecZeroEntries( tmp2 );
		VecAXPY( tmp1, 1.0, um0 );
		VecAXPY( tmp1, 0.5*ds, g0 );
		MatMult( A, tmp1, tmp2 );
		MatMult( Minv, tmp2, g1 );

		VecZeroEntries( tmp1 );
		VecZeroEntries( tmp2 );
		VecAXPY( tmp1, 1.0, um0 );
		VecAXPY( tmp1, 0.5*ds, g1 );
		MatMult( A, tmp1, tmp2 );
		MatMult( Minv, tmp2, g2 );

		VecZeroEntries( tmp1 );
		VecZeroEntries( tmp2 );
		VecAXPY( tmp1, 1.0, um0 );
		VecAXPY( tmp1, ds, g2 );
		MatMult( A, tmp1, tmp2 );
		MatMult( Minv, tmp2, g3 );

		VecZeroEntries( um1 );
		VecAXPY( um1, 1.0, um0 );
		VecAXPY( um1, ds/6.0, g0 );
		VecAXPY( um1, ds/3.0, g1 );
		VecAXPY( um1, ds/3.0, g2 );
		VecAXPY( um1, ds/6.0, g3 );

		VecCopy( um1, um0 );
	}

	VecToField( um1, phiTilde );

	VecDestroy( um0 );
	VecDestroy( um1 );
	VecDestroy( g0 );
	VecDestroy( g1 );
	VecDestroy( g2 );
	VecDestroy( g3 );
	VecDestroy( tmp1 );
	VecDestroy( tmp2 );

	delete psi;
}

void OIFS::SolveSecondOrder() {
	Field* 	field = new Field( "oifs-field", psi0->mesh, psi0->nDofs, psi0->bcs );
	Vec 	u0, u1, uHalf, uTmp1, uTmp2, um0, um1;
	int	cycle_i;
	double	ds = dt/nCycles;

	op->field = field;

	VecDuplicate( x, &u0 );
	VecDuplicate( x, &u1 );
	VecDuplicate( x, &uHalf );
	VecDuplicate( x, &uTmp1 );
	VecDuplicate( x, &uTmp2 );
	VecDuplicate( x, &um0 );
	VecDuplicate( x, &um1 );

	/* n - 1 time step */
	FieldToVec( phi1, u0 ); 

	for( cycle_i = 0; cycle_i < 2*nCycles; cycle_i++ ) {
		Extrapolate( field, -1.0 + (0.0 + 1.0*cycle_i)/nCycles );
		MatZeroEntries( A );
		matrix->Assemble();
		MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
	
		VecZeroEntries( uTmp1 );
		VecZeroEntries( uTmp2 );
		MatMult( A,    u0,    uTmp1 );
		MatMult( Minv, uTmp1, uTmp2 );
		VecZeroEntries( uHalf );
		VecAXPY( uHalf, 1.0,    u0    );
		VecAXPY( uHalf, 0.5*ds, uTmp2 );

		Extrapolate( field, -1.0 + (0.5 + 1.0*cycle_i)/nCycles );
		MatZeroEntries( A );
		matrix->Assemble();
		MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );

		VecZeroEntries( uTmp1 );
		VecZeroEntries( uTmp2 );
		MatMult( A,    uHalf, uTmp1 );
		MatMult( Minv, uTmp1, uTmp2 );
		VecZeroEntries( u1 );
		VecAXPY( u1, 1.0, u0    );
		VecAXPY( u1, ds,  uTmp2 );

		VecZeroEntries( u0 );
		VecCopy( u1, u0 );
	}

	VecZeroEntries( um1 );
	VecCopy( u1, um1 );

	/* n - 0 time step */
	FieldToVec( phi0, u0 );

	for( cycle_i = 0; cycle_i < nCycles; cycle_i++ ) {
		Extrapolate( field, 0.0 + (0.0 + cycle_i)/nCycles );
		MatZeroEntries( A );
		matrix->Assemble();
		MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );

		VecZeroEntries( uTmp1 );
		VecZeroEntries( uTmp2 );
		MatMult( A,    u0,    uTmp1 );
		MatMult( Minv, uTmp1, uTmp2 );
		VecZeroEntries( uHalf );
		VecAXPY( uHalf, 1.0,    u0    );
		VecAXPY( uHalf, 0.5*ds, uTmp2 );
	
		Extrapolate( field, 0.0 + (0.5 + cycle_i)/nCycles );
		MatZeroEntries( A );
		matrix->Assemble();
		MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
		MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );

		VecZeroEntries( uTmp1 );
		VecZeroEntries( uTmp2 );
		MatMult( A,    uHalf, uTmp1 );
		MatMult( Minv, uTmp1, uTmp2 );
		VecZeroEntries( u1 );
		VecAXPY( u1, 1.0, u0    );
		VecAXPY( u1, ds,  uTmp2 );

		VecZeroEntries( u0 );
		VecCopy( u1, u0 );
	}

	VecZeroEntries( um0 );
	VecCopy( u1, um0 );

	/* sum the contributions */
	VecZeroEntries( uTmp1 );
	VecAXPY( uTmp1, +4.0/3.0, um0 );
	VecAXPY( uTmp1, -1.0/3.0, um1 );
	//VecAXPY( uTmp1, +2.0, um0 );
	//VecAXPY( uTmp1, -0.5, um1 );
	VecToField( uTmp1, phiTilde );

	VecDestroy( u0 );
	VecDestroy( u1 );
	VecDestroy( uHalf );
	VecDestroy( uTmp1 );
	VecDestroy( uTmp2 );
	VecDestroy( um0 );
	VecDestroy( um1 );

	delete field;
}

void OIFS::Extrapolate( Field* field, double a ) {
	int node_i, dof_i;

	for( node_i = 0; node_i < psi0->mesh->nVertsTotal; node_i++ ) {
		for( dof_i = 0; dof_i < psi0->nDofs; dof_i++ ) {
			field->vals[node_i][dof_i] = (1.0 + a)*psi0->vals[node_i][dof_i] - a*psi1->vals[node_i][dof_i];
		}
	}
}

void InitVec( Vec* v, int size ) {
	VecCreate( MPI_COMM_WORLD, v );
	VecSetSizes( *v, size, PETSC_DETERMINE );
	VecSetType( *v, VECSEQ );
	VecSetOption( *v, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
	VecZeroEntries( *v );
}

void InitMat( Mat* A, int rSize, int cSize, int alloc ) {
	MatCreate( MPI_COMM_WORLD, A );
	MatSetSizes( *A, rSize, cSize, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( *A, MATSEQAIJ );
	MatSeqAIJSetPreallocation( *A, alloc, PETSC_NULL );
	MatZeroEntries( *A );
}

void FieldToVec( Field* field, Vec vec ) {
        int vec_i, field_i, dof_i;
        PetscScalar* array;
        int* map = field->bcs->fieldToVecMap;
        int nVertsTotal = field->mesh->nVertsTotal;

        VecGetArray( vec, &array );
        for( dof_i = 0; dof_i < field->nDofs; dof_i++ ) {
                for( field_i = 0; field_i < nVertsTotal; field_i++ ) {
                        vec_i = map[dof_i*nVertsTotal + field_i];
                        if( !field->bcs->IsBCNodeAndDof( field_i, dof_i ) ) {
                                array[vec_i] = field->vals[field_i][dof_i];
                        }
                }
        }
        VecRestoreArray( vec, &array );
}

void VecToField( Vec vec, Field* field ) {
        int vec_i, field_i, dof_i, map, vecSize;
        PetscScalar* array;
        int nVerts = field->mesh->nVertsTotal;
        int offset = field->bcs->vecOffset;

	VecGetSize( vec, &vecSize );
        VecGetArray( vec, &array );
        for( vec_i = 0; vec_i < vecSize; vec_i++ ) {
                map = field->bcs->vecToFieldMap[vec_i];
                field_i = map%nVerts;
                dof_i = map/nVerts;
                field->vals[field_i][dof_i] = array[offset + vec_i];
        }
        VecRestoreArray( vec, &array );

        /* update the periodic boundaries if required */
        if( field->mesh->periodic[0] || field->mesh->periodic[1] ) {
                field->PeriodicUpdate();
        }
}

void FieldCopy( Field* from, Field* to ) {
	int node_i, dof_i;

	for( node_i = 0; node_i < from->mesh->nVertsTotal; node_i++ ) {
		for( dof_i = 0; dof_i < from->nDofs; dof_i++ ) {
			if( !to->bcs->IsBCNodeAndDof( node_i, dof_i ) ) {
				to->vals[node_i][dof_i] = from->vals[node_i][dof_i];
			}
		}
	}
}
