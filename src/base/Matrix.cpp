#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"

using namespace std;
using std::string;

Matrix::Matrix( string _name, Mat _mat, Vector* _rowVector, Vector* _colVector, Vector* _rhsVector, Operator** _ops, int _nOps ) {
	name = _name;
	mat = _mat;
	rowVector = _rowVector;
	colVector = _colVector;
	rhsVector = _rhsVector;
	ops = _ops;
	nOps = _nOps;
}

Matrix::~Matrix() {
	for( int i = 0; i < nOps; i++ ) {
		delete ops[i];
	}
	delete[] ops;
}

void Matrix::AssembleElement( int el_i, double* M ) {
	for( int op_i = 0; op_i < nOps; op_i++ ) {
		ops[op_i]->AssembleElement( el_i, M );
	}
}

void Matrix::Assemble() {
	Mesh*		rowMesh 	= rowVector->field->mesh;
	Mesh*		colMesh 	= colVector->field->mesh;
	int*		rowNodes;
	int*		colNodes;
	int*		gRowNodes;
	int*		gColNodes;
	int		nRowNodes	= rowMesh->el->nNodes;
	int		nColNodes	= colMesh->el->nNodes;
	int		nRowDofs	= rowVector->field->nDofs;
	int		nColDofs	= colVector->field->nDofs;
	int		rowSize		= nRowNodes*nRowDofs;
	int		colSize		= nColNodes*nColDofs;
	int*		rowMap		= rowVector->field->bcs->fieldToVecMap;
	int*		colMap		= colVector->field->bcs->fieldToVecMap;
	double*		rhs;
	double* 	M;
	int		lRow, lCol, gRow, gCol;
	bool		addToVec;

	M = new double[rowSize*colSize];
	rhs = new double[rowSize];
	rowNodes = new int[rowSize];
	colNodes = new int[colSize];

	cout << "assembling " << name.c_str() << " matrix... ";

	for( int el_i = 0; el_i < rowMesh->nElsTotal; el_i++ ) {
		/* assemble the matrix */
		memset( M, 0, rowSize*colSize*sizeof(double) );
		AssembleElement( el_i, M );

		/* map the global node indices to matrix indices. note that indices that are entered in as '-1' (corresponding
		 * to dirichlet boundary condition nodes) will be ignored by petsc */
		gRowNodes = rowMesh->ElNodes( el_i );
		for( int node_i = 0; node_i < nRowNodes; node_i++ ) {
			for( int dof_i = 0; dof_i < nRowDofs; dof_i++ ) {
				lRow = node_i*nRowDofs + dof_i;
				gRow = dof_i*rowMesh->nVertsTotal + gRowNodes[node_i];
				rowNodes[lRow] = rowMap[gRow];
			}
		}

		gColNodes = colMesh->ElNodes( el_i );
		for( int node_i = 0; node_i < nColNodes; node_i++ ) {
			for( int dof_i = 0; dof_i < nColDofs; dof_i++ ) {
				lCol = node_i*nColDofs + dof_i;
				gCol = dof_i*colMesh->nVertsTotal + gColNodes[node_i];
				colNodes[lCol] = colMap[gCol];
			}
		}

		MatSetValues( mat, rowSize, rowNodes, colSize, colNodes, M, ADD_VALUES );

		/* copy any boundary conditions into the RHS vector (if specified) */
		if( !rhsVector ) {
			continue;
		}

		addToVec = false;
		memset( rhs, 0, rowSize*sizeof(double) );
		for( int node_j = 0; node_j < nColNodes; node_j++ ) {
			for( int dof_j = 0; dof_j < nColDofs; dof_j++ ) {
				lCol = node_j*nColDofs + dof_j;
				if( colNodes[lCol] == -1 ) {
					addToVec = true;
					for( int node_i = 0; node_i < nRowNodes; node_i++ ) {
						for( int dof_i = 0; dof_i < nRowDofs; dof_i++ ) {
							lRow = node_i*nRowDofs + dof_i;
							rhs[lRow] -= colVector->field->vals[gColNodes[node_j]][dof_j]*M[lRow*colSize + lCol];
						}
					}
				}
			}
		}
		if( addToVec ) {
			VecSetValues( rhsVector->vec, rowSize, rowNodes, rhs, ADD_VALUES );
		}
	}

	cout << "done.\n";

	delete[] M;
	delete[] rhs;
	delete[] rowNodes;
	delete[] colNodes;
}

void Matrix::Write( string filename, bool writeZeros ) {
        ofstream file;
        PetscInt row, nCols, col;
        const PetscInt* cols;
        const PetscScalar* vals;

        file.open( filename.c_str() );

        for( row = 0; row < (int)rowVector->vecSize; row++ ) {
                MatGetRow( mat, row, &nCols, &cols, &vals );
                for( col = 0; col < nCols; col++ ) {
			if( fabs(vals[col]) < 1.0e-8 ) {
				if( writeZeros ) {
					file << cols[col] << ": 0\t";
				}
			}
			else {
				file << cols[col] << ": " << vals[col] << "\t";
			}
                }
                file << "\n\n";
                MatRestoreRow( mat, row, &nCols, &cols, &vals );
        }

        file.close();
}

void Matrix::InverseDiagonal( Mat matInv ) {
	int 			size = 0, row, col, nCols, dof_i;
	const PetscInt* 	cols;
	const PetscScalar* 	vals;

	if( rowVector != colVector ) { 
		cerr << "ERROR: cannot take diagonal inverse of a matrix with different row and column variables" << endl;
		exit( 1 );
	}

	for( dof_i = 0; dof_i < rowVector->field->nDofs; dof_i++ ) {
		size += rowVector->field->mesh->nVertsTotal - rowVector->field->bcs->size[dof_i];
	}

	for( row = 0; row < size; row++ ) {
		MatGetRow( mat, row, &nCols, &cols, &vals );
		for( col = 0; col < nCols; col++ ) { if( cols[col] == row ) { break; } }
		MatSetValue( matInv, row, row, 1.0/vals[col], INSERT_VALUES );
		MatRestoreRow( mat, row, &nCols, &cols, &vals );
	}
	MatAssemblyBegin( matInv, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( matInv, MAT_FINAL_ASSEMBLY );
}
