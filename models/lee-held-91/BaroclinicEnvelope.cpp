#include <cmath>
#include <fstream>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "LinSys.h"
#include "TwoLayerQGEqn.h"

using namespace std;
using std::string;

void CalcTauAndChi( Field* psi1, Field* psi2, Field* tau, Field* chi ) {
	for( int node_i = 0; node_i < tau->mesh->nVertsTotal; node_i++ ) {
		tau->vals[node_i][0] = 0.5*(psi1->vals[node_i][0] - psi2->vals[node_i][0]);
		chi->vals[node_i][0] = 0.5*(psi1->vals[node_i][0] + psi2->vals[node_i][0]);
	}
}

void CalcPsi( Field* psi, Field* omega ) {
	PoissonEqn*	pe 	= new PoissonEqn( psi, omega, -1.0 );

	pe->Solve( false );
	delete pe;
}

double uTopIC( double* x ) { return exp( -0.1*(x[1]-10.5)*(x[1]-10.5) ); }
double omegaIC( double* x ) { return 0.2*(x[1]-10.5)*exp( -0.1*(x[1]-10.5)*(x[1]-10.5) ); }

int main( int argc, char** argv ) {
	char		tag[6]		= "petsc";
	int		nx[2]		= { 120, 20 };
	bool		periodic[2]	= { true, false };
	double		min[2]		= { 0, 0 };
	double		max[2]		= { 40.0*M_PI, 21.0 };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", 4, min, max, periodic );
	BCs*		bcs		= new BCs( true, true, false, false, mesh, 0 );
	Field*		psi1		= new Field( "psi-1-curr", mesh, 1, bcs );
	Field*		psi2		= new Field( "psi-2-curr", mesh, 1, bcs );
	Field*		omega1		= new Field( "omega-1-curr", mesh, 1, bcs );
	Field*		omega2		= new Field( "omega-2-curr", mesh, 1, bcs );
	Field*		psi1Prev	= new Field( "psi-1-prev", mesh, 1, NULL );
	Field*		psi2Prev	= new Field( "psi-2-prev", mesh, 1, NULL );
	Field*		omega1Prev	= new Field( "omega-1-prev", mesh, 1, NULL );
	Field*		omega2Prev	= new Field( "omega-2-prev", mesh, 1, NULL );
	Field*		tau		= new Field( "tau", mesh, 1, NULL );
	Field*		chi		= new Field( "chi", mesh, 1, NULL );
	Field*		tau_e		= new Field( "tau-e", mesh, 1, bcs );
	Field*		fields[8];
	double		dt		= 0.01;
	int		timeStep	= atoi( argv[1] );
	double		time		= timeStep*dt;
	int		dumpEvery	= 40;
	int		nTimeSteps	= 40000;
	TwoLayerQGEqn*	qg;

	PetscInitialize( &argc, &argv, (char)0, tag );

	if( timeStep > 0 ) {
		psi1->Read( timeStep );
		psi2->Read( timeStep );
		omega1->Read( timeStep );
		omega2->Read( timeStep );
	}
	else {
		omega1->SetICFunc( 0, omegaIC );
		CalcPsi( psi1, omega1 );
	}
	fields[0] = psi1;
	fields[1] = psi2;
	fields[2] = omega1;
	fields[3] = omega2;
	fields[4] = tau;
	fields[5] = chi;

	mesh->Save();

	CalcPsi( tau_e, omega1 );

	qg = new TwoLayerQGEqn( psi1, omega1, psi2, omega2, 0.5, 0.5, 0.30, 0.05, 0.033333, 0.0001, tau_e, dt );

	if( timeStep == 0 ) {
		CalcTauAndChi( psi1, psi2, tau, chi );
		fields[6] = qg->GenVel( psi1, "vel-1" );
		fields[7] = qg->GenVel( psi2, "vel-2" );
		WriteXDMFHeader( timeStep );
		WriteXDMF( fields, 8, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
		delete fields[6];
		delete fields[7];
	}

	psi1Prev->Copy( psi1 );
	psi2Prev->Copy( psi2 );
	omega1Prev->Copy( omega1 );
	omega2Prev->Copy( omega2 );

	time += dt;
	timeStep++;
	qg->Solve( NULL, NULL, NULL, NULL );

	CalcTauAndChi( psi1, psi2, tau, chi );
	fields[6] = qg->GenVel( psi1, "vel-1" );
	fields[7] = qg->GenVel( psi2, "vel-2" );
	psi1Prev->Copy( psi1 );
	psi2Prev->Copy( psi2 );
	omega1Prev->Copy( omega1 );
	omega2Prev->Copy( omega2 );
	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 8, timeStep, time, dt );
	WriteXDMFFooter( timeStep );
	delete fields[6];
	delete fields[7];

	for( timeStep++; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		qg->Solve( psi1Prev, omega1Prev, psi2Prev, omega2Prev );
		if( timeStep%dumpEvery == 0 ) {
			fields[6] = qg->GenVel( psi1, "vel-1" );
			fields[7] = qg->GenVel( psi2, "vel-2" );
			CalcTauAndChi( psi1, psi2, tau, chi );
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 8, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
			delete fields[6];
			delete fields[7];
		}
	}

	delete qg;
	delete mesh;
	delete bcs;
	delete psi1;
	delete psi2;
	delete omega1;
	delete omega2;
	delete psi1Prev;
	delete psi2Prev;
	delete omega1Prev;
	delete omega2Prev;
	delete tau;
	delete chi;
	delete tau_e;

	PetscFinalize();

	return EXIT_SUCCESS;
}
