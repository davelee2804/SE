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
#include "PoissonEqn.h"

using namespace std;
using std::string;

double cosineBC( double* coord ) { return                cos( M_PI*coord[0] )*cos( M_PI*coord[1] ); }
double cosineIC( double* coord ) { return -2.0*M_PI*M_PI*cos( M_PI*coord[0] )*cos( M_PI*coord[1] ); }

double IntegrateError( Field* psi ) {
	int el_i, pt_i;
	double errorSq = 0.0, analyticSq = 0.0, numeric, analytic, gCoord[2];
	double* lCoord, weight, detJac;

	for( el_i = 0; el_i < psi->mesh->nElsTotal; el_i++ ) {
		for( pt_i = 0; pt_i < psi->mesh->el->nPoints; pt_i++ ) {
			lCoord = psi->mesh->el->quadPts[pt_i]->coord;
			weight = psi->mesh->el->quadPts[pt_i]->weight;
			detJac = psi->mesh->DetJac( el_i, pt_i );
			psi->InterpLocal( el_i, lCoord, &numeric );
			psi->mesh->LocalToGlobal( lCoord, el_i, gCoord );
			analytic = cosineBC( gCoord );
			errorSq += detJac*weight*( analytic - numeric )*( analytic - numeric );
			analyticSq += detJac*weight*analytic*analytic;
		}
	}

	return sqrt( errorSq / analyticSq );
}

int main( int argc, char** argv ) {
	int		N		= atoi( argv[1] );
	int 		nx[2] 		= { 4, 4 };
	double 		min[2] 		= { -1.0, -1.0 };
	double 		max[2] 		= { +1.0, +1.0 };
	bool 		periodic[2] 	= { false, false };
	Mesh* 		mesh		= new Mesh( "mesh", nx, "legendre", N, min, max, periodic );
	BCs* 		bcs		= new BCs( true, true, true, true, mesh, 0 );
	Field* 		psi		= new Field( "psi", mesh, 1, bcs );
	Field* 		omega		= new Field( "omega", mesh, 1, NULL );
	Field*		analytic	= new Field( "analytic", mesh, 1, NULL );
	PoissonEqn* 	poissonEqn;
	Field* 		fields[3];
	char 		tag[] 		= "petsc";

	PetscInitialize( &argc, &argv, (char*)0, tag );

	/* setup */
	psi->SetBCFunc( "bottom", 0, cosineBC );
	psi->SetBCFunc( "top",    0, cosineBC );
	psi->SetBCFunc( "left",   0, cosineBC );
	psi->SetBCFunc( "right",  0, cosineBC );
	omega->SetICFunc( 0, cosineIC );
	analytic->SetICFunc( 0, cosineBC );

	poissonEqn = new PoissonEqn( psi, omega, -1.0 );

	/* solve */
	poissonEqn->Solve( false );

	/* write results */
	mesh->Save();
	fields[0] = psi;
	fields[1] = analytic;
	fields[2] = omega;
	WriteXDMFHeader( 1 );
	WriteXDMF( fields, 3, 1, 0.0, 1 );
	WriteXDMFFooter( 1 );

	cout << "polynomial order: " << mesh->N << " error: " << IntegrateError( psi ) << endl;

	/* cleanup */
	delete poissonEqn;
	delete analytic;
	delete omega;
	delete psi;
	delete bcs;
	delete mesh;	

	PetscFinalize();

	return 0;
}
