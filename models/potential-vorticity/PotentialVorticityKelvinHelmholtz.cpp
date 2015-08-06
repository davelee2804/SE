#include <iostream>
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
#include "Utils.h"
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "Advector.h"
#include "PoissonEqn.h"

using namespace std;
using std::string;

double omegaIC( Coord coord ) {
        double mid = 0.5;
        double a = 0.08454459;
        double U_0 = 1.0;

        return (fabs(coord[1] - mid) <= a) ? U_0/a : 0.0;
}

void GetVel( Field* psi, Field* vel ) {
        int node_i;
        double **gradPsi;

	gradPsi = new double*[1];
	gradPsi[0] = new double[2];

        for( node_i = 0; node_i < psi->mesh->nVertsTotal; node_i++ ) {
                psi->InterpDerivsGlobal( psi->mesh->verts[node_i], gradPsi );
                vel->vals[node_i][0] = -gradPsi[0][1];
                vel->vals[node_i][1] = +gradPsi[0][0];
        }
	delete[] gradPsi[0];
	delete[] gradPsi;
}

void UpdateOmega( Field* omega, Field* omegaSL1, Field* omegaSL2 ) {
	int node_i;

	for( node_i = 0; node_i < omega->mesh->nVertsTotal; node_i++ ) {
		if( !omega->bcs->IsBCNode( node_i ) ) {
			omega->vals[node_i][0] = (4.0/3.0)*omegaSL1->vals[node_i][0] - (1.0/3.0)*omegaSL2->vals[node_i][0];
			//omega->vals[node_i][0] = omegaSL1->vals[node_i][0];
		}
	}

	omega->PeriodicUpdate();
}

int main( int argc, char** argv ) {
	char 		tag[]		= "petsc";
	int		N		= 5;
	int		nx[2]		= { 80, 20 };
	double		min[2]		= { 0.0, 0.0 };
	double		max[2]		= { 4.0, 1.0 };
	bool		periodic[2]	= { true, false };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", N, min, max, periodic );
	BCs*		bcs		= new BCs( true, true, false, false, mesh, 0 );
	Field*		omega		= new Field( "omega", mesh, 1, bcs );
	Field*		psi		= new Field( "psi", mesh, 1, bcs );
	Field*		vel		= new Field( "vel", mesh, 2, NULL );
	Field*		prevOmega	= new Field( "prevOmega", mesh, 1, bcs );
	Field*		prevVel		= new Field( "prevVel", mesh, 2, NULL );
	Field*		fields[3];
	double		dt		= 1000.0*pow( ((max[0] - min[0])/nx[0])/N, ((double)(N + 1))/(2 + 1));
	double		time		= 0.0;
	int		nTimeSteps	= 10000;
	int		timeStep_i;
	int		dumpEvery	= 25;
	Advector*	advector;
	PoissonEqn*	poisson;

	PetscInitialize( &argc, &argv, (char*)0, tag );

	/* initial and boundary conditions */
	omega->SetICFunc( 0, omegaIC );
	psi->SetBCConst( "bottom", 0, -0.9577277 );
	psi->SetBCConst( "top", 0, -0.9577277 );

	/* initial solve for psi */
	poisson = new PoissonEqn( psi, omega );
	poisson->Solve();
	psi->PeriodicUpdate();
	delete poisson;
	GetVel( psi, vel );

	/* write the initial fields to file */
	WriteXDMFHeader( 0 );
	fields[0] = omega;
	fields[1] = psi;
	fields[2] = vel;
	WriteXDMF( fields, 3, 0, time, 0.0 );
	WriteXDMFFooter( 0 );

	prevOmega->Copy( omega );
	prevVel->Copy( vel );

	for( timeStep_i = 1; timeStep_i <= nTimeSteps; timeStep_i++ ) {
		time += dt;
		cout << "time step: " << timeStep_i << ", time: " << time << ", dt: " << dt << endl;

		/* s-L advect omega */
		advector = new Advector( omega, vel, prevOmega, prevVel );
		advector->Advect( dt );
		prevOmega->Copy( omega );
		UpdateOmega( omega, advector->fieldSL, advector->fieldSLMinusOne );
		delete advector;

		/* solve the eliptic equation */
		poisson = new PoissonEqn( psi, omega );
		poisson->Solve();
		psi->PeriodicUpdate();
		prevVel->Copy( vel );
		GetVel( psi, vel );
		delete poisson;

		if( timeStep_i%dumpEvery == 0 ) {
			WriteXDMFHeader( timeStep_i );
			WriteXDMF( fields, 3, timeStep_i, time, dt );
			WriteXDMFFooter( timeStep_i );
		}
	}

	PetscFinalize();

	delete prevVel;
	delete prevOmega;
	delete vel;
	delete psi;
	delete omega;
	delete bcs;
	delete mesh;

	return 0;
}
