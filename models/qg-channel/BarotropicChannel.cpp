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
#include "QuasiGeostrophicEqn.h"

using namespace std;
using std::string;

#define BETA 0.0
#define F 0.0

#define U0 1.0
#define HALF_WIDTH 1.0

//double stream_func_ic( double* x ) { return (x[1]*x[1] - M_PI*M_PI)*pow( x[1], 3 ); }
//double vorticity_ic( double* x ) { return 20.0*pow( x[1], 3 ) - 6.0*M_PI*M_PI*x[1]; }
//double stream_func_ic( double* x ) { return (x[1]*x[1] - M_PI*M_PI)*pow( x[1], 4 ); }
//double vorticity_ic( double* x ) { return 30.0*pow( x[1], 4 ) - 12.0*M_PI*M_PI*pow( x[1], 2 ); }
//double stream_func_ic( double* x ) { return 1.0 - ( 32*pow( x[1]/M_PI, 6 ) - 48*pow( x[1]/M_PI, 4 ) + 18*pow( x[1]/M_PI, 2 ) - 1 ); }
//double vorticity_ic( double* x ) { return (-1.0/(M_PI*M_PI))*( 32*6*5*pow( x[1]/M_PI, 4 ) - 48*4*3*pow( x[1]/M_PI, 2 ) + 18*2*1 ); }
//double stream_func_ic( double* x ) { return 1.0 - ( 8*pow( x[1]/M_PI, 4 ) - 8*pow( x[1]/M_PI, 2 ) + 1 ); }
//double vorticity_ic( double* x ) { return (-1.0/(M_PI*M_PI))*( 96*pow( x[1]/M_PI, 2 ) - 16 ); }

/*double stream_func_ic( double* x ) { 
	if( x[1] > +HALF_WIDTH )      { return -U0*x[1] + 0.5*U0*HALF_WIDTH; }
	else if( x[1] < -HALF_WIDTH ) { return +U0*x[1] + 0.5*U0*HALF_WIDTH; }
	else                          { return -0.5*U0*x[1]*x[1]/HALF_WIDTH; }
}
double vorticity_ic( double* x ) {
	if( fabs( x[1] ) > HALF_WIDTH ) { return 0.0; }
	else                            { return -U0/HALF_WIDTH; }
}*/

/*double vorticity_ic( double* x ) { return tanh( 10.0*x[1] ); }
void CalcStreamFunc( Field* phi, Field* omega ) {
	PoissonEqn* 	poisson		= new PoissonEqn( phi, omega, -1.0 );

	poisson->Solve( false );

	delete poisson;
}*/
#define A0 8
double vorticity_ic( double* x )   { return -A0/( cosh( A0*x[1] )*cosh( A0*x[1] ) ); }
double stream_func_ic( double* x ) { return -log( cosh( A0*x[1] ) )/A0; }

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 48, 12 };
	bool			periodic[2]	= { true, false };
	double			min[2]		= { -4*M_PI, -M_PI };
	double			max[2]		= { +4*M_PI, +M_PI };
	Mesh*			mesh		= new Mesh( "mesh", nx, "legendre", 8, min, max, periodic );
	BCs*			bcs		= new BCs( true, true, false, false, mesh, 0 );
	BCs*			none		= new BCs( false, false, false, false, mesh, 0 );
	Field*			phi		= new Field( "phi", mesh, 1, bcs );
	Field*			omega		= new Field( "omega", mesh, 1, none );
	Field*			phiPrev		= new Field( "phi-prev", mesh, 1, bcs );
	Field*			omegaPrev	= new Field( "omega-prev", mesh, 1, NULL );
	Field*			fields[3];
	double			dt		= 0.5*(8*M_PI/nx[0])/U0;
	int			timeStep	= atoi( argv[1] );
	double			time		= timeStep*dt;
	int			dumpEvery	= 10;
	int			nTimeSteps	= 8000;
	QuasiGeostrophicEqn*	qg;

	PetscInitialize( &argc, &argv, (char)0, tag );

	mesh->Save();

	if( timeStep > 0 ) {
		phi->Read( timeStep );
		omega->Read( timeStep );
	}
	else {
		phi->SetICFunc( 0, stream_func_ic );
		omega->SetICFunc( 0, vorticity_ic );
		//CalcStreamFunc( phi, omega );
	}
	qg = new QuasiGeostrophicEqn( phi, omega, F, BETA, 0.0004, 0.0, dt, NULL );
	fields[0] = phi;
	fields[1] = omega;
	fields[2] = qg->GenVel( phi );

	if( timeStep == 0 ) {
		WriteXDMFHeader( timeStep );
		WriteXDMF( fields, 3, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
	}

	/* first time step */
	time += dt;
	timeStep++;

	phiPrev->Copy( phi );
	omegaPrev->Copy( omega );

	qg->Solve( NULL, NULL );

	delete fields[2];
	fields[2] = qg->GenVel( phi );
	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 3, timeStep, time, dt );
	WriteXDMFFooter( timeStep );

	for( timeStep++; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		qg->Solve( phiPrev, omegaPrev );
		if( timeStep%dumpEvery == 0 ) {
			delete fields[2];
			fields[2] = qg->GenVel( phi );
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 3, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
		}
	}

	delete qg;
	delete mesh;
	delete bcs;
	delete phi;
	delete omega;
	delete phiPrev;
	delete omegaPrev;
	delete fields[2];

	PetscFinalize();

	return EXIT_SUCCESS;
}
