#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "QuadPoint.h"
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
#include "PVEqn.h"

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

int main( int argc, char** argv ) {
	char		tag[]		= "petsc";
	int		N		= 7;
	int		nx[2]		= { 32, 8 };
	double		min[2]		= { 0.0, 0.0 };
	double		max[2]		= { 4.0, 1.0 };
	bool		periodic[2]	= { true, false };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", N, min, max, periodic );
	BCs*		omegaBCs	= new BCs( true, true, false, false, mesh, 0 );
	int		omegaSize	= mesh->nVertsTotal - omegaBCs->size[0];
	BCs*		psiBCs		= new BCs( true, true, false, false, mesh, omegaSize );
	BCs*		tempPsiBCs	= new BCs( true, true, false, false, mesh, 0 );
	Field*		omega		= new Field( "omega", mesh, 1, omegaBCs );
	Field*		psi		= new Field( "psi", mesh, 1, psiBCs );
	Field*		tempPsi		= new Field( "tempPsi", mesh, 1, tempPsiBCs );
	Field*		prevOmega	= new Field( "prevOmega", mesh, 1, omegaBCs );
	Field*		tempOmega	= new Field( "tempOmega", mesh, 1, omegaBCs );
	Field*		velocity	= new Field( "velocity", mesh, 2, NULL );
	Field*		prevVel		= new Field( "prevVel", mesh, 2, NULL );
	PoissonEqn*	poisson;
	PVEqn*		pv;
	double		time		= 0.0;
	double          dt              = 0.1*(max[0] - min[0])/nx[0];//pow( ((max[0] - min[0])/nx[0])/N, ((double)(N + 1))/(2 + 1));
	int		timeStep_i;
	int		nTimeSteps	= 10000;
	int 		dumpEvery	= 1;
	Field*		fields[3];

	PetscInitialize( &argc, &argv, (char*)0, tag );

	/* set the initial and boundary conditions */
	omega->SetICFunc( 0, omegaIC );
	tempPsi->SetBCConst( "bottom", 0, -0.9577277 );
	tempPsi->SetBCConst( "top", 0, -0.9577277 );

	/* initial guess for psi */
	poisson = new PoissonEqn( tempPsi, omega );
	poisson->Solve();
	tempPsi->PeriodicUpdate();
	delete poisson;
	GetVel( tempPsi, velocity );

	fields[0] = omega;
	fields[1] = tempPsi;
	fields[2] = velocity;

	/* write the initial values */
	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 3, 0, time, 0.0 );
	WriteXDMFFooter( 0 );

	/* first time step, use single velocity & omega for advection only */
	psi->Copy( tempPsi );
	prevOmega->Copy( omega );
	time += dt;
	cout << "time step: " << 1 << ", time: " << time << ", dt: " << dt << endl;
	pv = new PVEqn( omega, psi, velocity, NULL, NULL );
	pv->Solve( dt );
	delete pv;
	omega->PeriodicUpdate();
	psi->PeriodicUpdate();
	prevVel->Copy( velocity );
	GetVel( psi, velocity );

	fields[1] = psi;

	if( dumpEvery == 1 ) {
		WriteXDMFHeader( 1 );
		WriteXDMF( fields, 3, 1, time, dt );
		WriteXDMFFooter( 1 );
	}

	for( timeStep_i = 2; timeStep_i <= nTimeSteps; timeStep_i++ ) {
		time += dt;
		cout << "time step: " << timeStep_i << ", time: " << time << ", dt: " << dt << endl;
		tempOmega->Copy( omega );
		pv = new PVEqn( omega, psi, velocity, prevOmega, prevVel );
		pv->Solve( dt );
		delete pv;
		omega->PeriodicUpdate();
		psi->PeriodicUpdate();
		prevOmega->Copy( tempOmega );
		prevVel->Copy( velocity );
		GetVel( psi, velocity );

		if( timeStep_i%dumpEvery == 0 ) {
			WriteXDMFHeader( timeStep_i );
			WriteXDMF( fields, 3, timeStep_i, time, dt );
			WriteXDMFFooter( timeStep_i );
		}
	}

	delete prevVel;
	delete velocity;
	delete tempOmega;
	delete prevOmega;
	delete tempPsi;
	delete psi;
	delete omega;
	delete tempPsiBCs;
	delete psiBCs;
	delete omegaBCs;
	delete mesh;

	PetscFinalize();

	return 0;
}
