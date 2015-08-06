#include <string>
#include <cmath>
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
#include "Utils.h"
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "HelmholtzEqn.h"

using namespace std;
using std::string;

#define K0 0.5*M_PI

double f( double* x ) { return cos( K0*x[0] )*cos( K0*x[1] ); }
double b( double* x ) { return (1.0 - 2.0*K0*K0)*cos( K0*x[0] )*cos( K0*x[1] ); }

double IntegrateError( Field* psi, int dof, BCFunc* func ) {
	double errorSq = 0.0, analyticSq = 0.0, numeric[2], analytic, gCoord[2];
	double* lCoord, weight, detJac;

	for( int el_i = 0; el_i < psi->mesh->nElsTotal; el_i++ ) {
		for( int pt_i = 0; pt_i < psi->mesh->el->nPoints; pt_i++ ) {
			lCoord = psi->mesh->el->quadPts[pt_i]->coord;
			weight = psi->mesh->el->quadPts[pt_i]->weight;
			detJac = psi->mesh->DetJac( el_i, pt_i );
			psi->InterpLocal( el_i, lCoord, numeric );
			psi->mesh->LocalToGlobal( lCoord, el_i, gCoord );
			analytic    = func( gCoord );
			errorSq    += detJac*weight*( analytic - numeric[dof] )*( analytic - numeric[dof] );
			analyticSq += detJac*weight*analytic*analytic;
		}
	}

	return sqrt( errorSq / analyticSq );
}

int main( int argc, char** argv ) {
	int		N		= 8;
	int 		nx[2] 		= { 4, 4 };
	double 		min[2] 		= { -1.0, -1.0 };
	double 		max[2] 		= { +1.0, +1.0 };
	bool 		periodic[2] 	= { false, false };
	Mesh* 		mesh		= new Mesh( "mesh", nx, "chebyshev", N, min, max, periodic );
	bool		noslip[2]	= { true, true };
	BCs* 		bcs		= new BCs( noslip, noslip, noslip, noslip, mesh, 0 );
	Field* 		psi		= new Field( "psi", mesh, 1, bcs );
	Field*		analytic	= new Field( "analytic", mesh, 1, NULL );
	Field*		omega		= new Field( "omega", mesh, 1, NULL );
	HelmholtzEqn* 	helmholtzEqn;
	Field* 		fields[3];
	char 		tag[] 		= "petsc";

	PetscInitialize( &argc, &argv, (char*)0, tag );

	/* setup */
	analytic->SetICFunc( 0, f );
	omega->SetICFunc( 0, b );

	helmholtzEqn = new HelmholtzEqn( psi, omega, 1.0, -1.0, 1.0, 0 );

	/* write results */
	fields[0] = psi;
	fields[1] = analytic;
	fields[2] = omega;
	mesh->Save();
	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 3, 0, 0.0, 1 );
	WriteXDMFFooter( 0 );

	/* solve */
	helmholtzEqn->Solve( "helm_", false );

	WriteXDMFHeader( 1 );
	WriteXDMF( fields, 3, 1, 0.0, 1 );
	WriteXDMFFooter( 1 );

	cout << "polynomial order: " << mesh->N << "\terror: " << IntegrateError( psi, 0, f ) << endl;

	/* cleanup */
	delete helmholtzEqn;
	delete analytic;
	delete omega;
	delete psi;
	delete bcs;
	delete mesh;	

	PetscFinalize();

	return 0;
}
