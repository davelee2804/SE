#include <string>
#include <iostream>
#include <fstream>
#include <petsc.h>
#include <petscvec.h>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "RHSOp.h"
#include "Vector.h"

using namespace std;
using std::string;

Vector::Vector( string _name, Field* _field, Vec _vec, RHSOp** _ops, int _nOps ) {
	name = _name;
	field = _field;
	vec = _vec;
	vecSize = field->mesh->nVertsTotal*field->nDofs;
	if( field->bcs ) {
		for( int dof_i = 0; dof_i < field->nDofs; dof_i++ ) {
			vecSize -= field->bcs->size[dof_i];
		}
	}
	ops = _ops;
	nOps = _nOps;
}

Vector::~Vector() {
	for( int i = 0; i < nOps; i++ ) {
		delete ops[i];
	}
	delete[] ops;
}

void Vector::UpdateField() {
	int vec_i, field_i, dof_i, map;
	PetscScalar* array;
	int nVerts = field->mesh->nVertsTotal;
	int offset = field->bcs->vecOffset;

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

void Vector::UpdateVec() {
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

void Vector::AssembleElement( int el_i, double* rhs ) {
	for( int op_i = 0; op_i < nOps; op_i++ ) {
		ops[op_i]->AssembleElement( el_i, rhs );
	}
}

void Vector::Assemble() {
	int		nElNodes  	= field->mesh->el->nNodes;
	int		elVecSize 	= nElNodes*field->nDofs;
	int		nVertsTotal	= field->mesh->nVertsTotal;
	double*		rhs;
	int*		vecNodes;
	int*		gNodes;
	int*		map		= field->bcs->fieldToVecMap;

	rhs = new double[elVecSize];
	vecNodes = new int[elVecSize];

	cout << "assembling " << name.c_str() << " vector... ";

	for( int el_i = 0; el_i < field->mesh->nElsTotal; el_i++ ) {
		memset( rhs, 0, elVecSize*sizeof(double) );
		AssembleElement( el_i, rhs );

		gNodes = field->mesh->ElNodes( el_i );
		for( int node_i = 0; node_i < nElNodes; node_i++ ) {
			for( int dof_i = 0; dof_i < field->nDofs; dof_i++ ) {
				vecNodes[node_i*field->nDofs + dof_i] = map[dof_i*nVertsTotal + gNodes[node_i]];
			}
		}

		VecSetValues( vec, elVecSize, vecNodes, rhs, ADD_VALUES );
	}

	cout << "done.\n";

	delete[] rhs;
	delete[] vecNodes;
}

void Vector::CopyTo( Vec v ) {
	PetscScalar *from, *to;
	int offset = field->bcs->vecOffset;

	VecGetArray( v, &to );
	VecGetArray( vec, &from );

	for( int i = offset; i < offset + vecSize; i++ ) {
		to[i] = from[i];
	}

	VecRestoreArray( v, &to );
	VecRestoreArray( vec, &from );
}

void Vector::CopyFrom( Vec v ) {
	PetscScalar *from, *to;
	int offset = field->bcs->vecOffset;

	VecGetArray( v, &from );
	VecGetArray( vec, &to );

	for( int i = offset; i < offset + vecSize; i++ ) {
		to[i] = from[i];
	}

	VecRestoreArray( v, &from );
	VecRestoreArray( vec, &to );
}

void Vector::Write( string filename ) {
	ofstream file;
	PetscScalar* array;
	double value;
	int offset = field->bcs->vecOffset;

	VecGetArray( vec, &array );

	file.open( filename.c_str() );

	for( int i = offset; i < offset + vecSize; i++ ) {
		value = (fabs(array[i]) < 1.0e-8) ? 0.000000 : array[i];
		file << value << endl;
	}

	file.close();

	VecRestoreArray( vec, &array );
}

double Vector::Norm() {
	int offset = field->bcs->vecOffset;
	double normSq = 0.0;
	PetscScalar* array;
	
	VecGetArray( vec, &array );

	for( int i = offset; i < offset + vecSize; i++ ) {
		normSq += array[i]*array[i];
	}

	VecRestoreArray( vec, &array );

	return sqrt( normSq );
}
