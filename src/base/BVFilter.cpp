#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Legendre.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "BVFilter.h"

using namespace std;
using std::string;

BVFilter::BVFilter( Field* _nodalField, double _nu, int _p, int _s ) {
	Legendre* l;

	nodalField = _nodalField;
	nu = _nu;
	p = _p;
	s = _s;

	modalField = new Field( "modalField", nodalField->mesh, nodalField->nDofs, nodalField->bcs );

	l = (Legendre*)nodalField->mesh->el;

	Binv = new double[(nodalField->mesh->el->nNodes)*(nodalField->mesh->el->nNodes)];
	w_k  = new double[nodalField->mesh->el->N+1];

	for( int k = 0; k <= nodalField->mesh->el->N; k++ ) { 
		w_k[k] = Sigma( k ); 
		cout << k << ":\t" << w_k[k] << endl;
	}

	B = l->ModalToNodalTransformMatrix();
	Inv( B, Binv, nodalField->mesh->el->nNodes );
}

BVFilter::~BVFilter() {
	delete modalField;
	delete[] B;
	delete[] Binv;
	delete[] w_k;
}

double BVFilter::Sigma( int k ) {
	int	N	= nodalField->mesh->el->N;
        double  theta, ki, omega, sigma;

        if( k < s ) { return 1.0; }

        theta = ((double)(k - s))/(N - s);
        omega = fabs( theta ) - 0.5;
        ki    = ( 2*(k - s) == (N - s) ) ? 1.0 : sqrt( -log(1.0 - 4.0*omega*omega)/(4.0*omega*omega) );
        sigma = 0.5*erfc( 2.0*sqrt( (double)p )*ki*omega );

        return sigma;
}

void BVFilter::Apply() {
	Mesh*		mesh		= nodalField->mesh;
	int		nNodes		= mesh->el->nNodes;
	int*		elNodes;
	double*		nodalElVec	= new double[nNodes];
	double*		modalElVec	= new double[nNodes];
	double*		filtElVec	= new double[nNodes];

	for( int el_i = 0; el_i < mesh->nElsTotal; el_i++ ) {
		elNodes = mesh->ElNodes( el_i );
		for( int dof_i = 0; dof_i < nodalField->nDofs; dof_i++ ) {
			for( int node_i = 0; node_i < nNodes; node_i++ ) {
				nodalElVec[node_i] = nodalField->vals[elNodes[node_i]][dof_i];
			}
			AXEB( Binv, nodalElVec, modalElVec, nNodes );
			for( int node_i = 0; node_i < nNodes; node_i++ ) {
				modalElVec[node_i] *= w_k[node_i%(mesh->el->N+1)]*w_k[node_i/(mesh->el->N+1)];
			}
			AXEB( B, modalElVec, filtElVec, nNodes );
			for( int node_i = 0; node_i < nNodes; node_i++ ) {
				nodalField->vals[elNodes[node_i]][dof_i] = filtElVec[node_i];
			}
		}
	}

	delete[] nodalElVec;
	delete[] modalElVec;
	delete[] filtElVec;
}
