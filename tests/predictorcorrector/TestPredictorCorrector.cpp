#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SW.h"
#include "2LRL.h"

using namespace std;
using std::string;

#define A 0.771
#define B 0.394

double u0( double* x ) { return A*B*B*0.25*( 6.0*x[1]*x[1] - 9.0 )*exp( -0.5*x[1]*x[1] )/( cosh( B*x[0] )*cosh( B*x[0] ) ); }
double v0( double* x ) { return -4.0*A*B*B*B*x[1]*tanh( B*x[0] )  *exp( -0.5*x[1]*x[1] )/( cosh( B*x[0] )*cosh( B*x[0] ) ); }
double h0( double* x ) { return A*B*B*0.25*( 6.0*x[1]*x[1] + 3.0 )*exp( -0.5*x[1]*x[1] )/( cosh( B*x[0] )*cosh( B*x[0] ) ); }

void SolvePC( PredictorCorrector* pc ) {
	Vec 	f, x, f_utn, f_ein, f_dhun, f_dhup, f_utnln, f_utnlp, f_utlp;
	Vector*	sol;
	Field*	velTopTemp = new Field( "vel-top-temp", pc->velTop->mesh, 2, NULL );
	double	alpha      = 2.0/3.0;

	/* vector steup */
	CreateVector( pc->velTop, &f_utn );		/* previous top layer velocity */
	CreateVector( pc->etaInt, &f_ein );		/* previous height interface perturbation */
	CreateVector( pc->etaInt, &f_dhun );	/* mass eqn nonlinear predictor */
	CreateVector( pc->etaInt, &f_dhup );	/* mass eqn nonlinear currector */
	CreateVector( pc->velTop, &f_utnln );	/* top layer non linear predictor */
	CreateVector( pc->velTop, &f_utnlp );	/* top layer non linear corrector */
	CreateVector( pc->velTop, &f_utlp );	/* top layer linear corrector */

	pc->advTopP->Advect( pc->dt );
	velTopTemp->Copy( pc->velTop );

	/* top layer momentum eqn predictor step */
	CreateVector( pc->velTop, &f );
	CreateVector( pc->velTop, &x );
	sol = new Vector( "vel-top-pred-sol", pc->velTopP, x, NULL, 0 );
	pc->AssemblePrevTimeStep( pc->advTopP->fieldSL, f_utn );
	pc->AssembleVelNonLinearRHS( pc->velTop, pc->etaInt, pc->tp, f_utnln );
	VecAXPY( f, 1.0, f_utn );
	VecAXPY( f, 1.0*alpha, f_utnln );
	VecAXPY( f, 1.0, pc->b_ut );
	pc->SolveLinSys( "tlp_", pc->velTop, pc->A1N, f, x );
	sol->UpdateField();
	VecDestroy( f );
	VecDestroy( x );
	delete sol;

	/* mass eqn predictor step */
	CreateVector( pc->etaInt, &f );
	CreateVector( pc->etaInt, &x );
	sol = new Vector( "eta-int-pred-sol", pc->etaIntP, x, NULL, 0 );
	pc->AssemblePrevTimeStep( pc->etaInt, f_ein );
	pc->AssembleMassEqnRHS( pc->velBot, pc->etaInt, f_dhun );
	VecAXPY( f, 1.0, f_ein );
	VecAXPY( f, 1.0, f_dhun );
	VecAXPY( f, 1.0, pc->b_ei );
	pc->SolveLinSys( "ihp_", pc->etaInt, pc->Ae, f, x );
	sol->UpdateField();
	VecDestroy( f );
	VecDestroy( x );
	delete sol;

	/* top layer momentum eqn corrector step */
	CreateVector( pc->velTop, &f );
	CreateVector( pc->velTop, &x );
	sol = new Vector( "vel-top-corr-sol", pc->velTop, x, NULL, 0 );
	pc->AssembleVelNonLinearRHS( pc->velTopP, pc->etaIntP, pc->tp, f_utnlp );
	pc->AssembleVelLinearRHS( pc->velTopP, pc->tp, f_utlp );
	VecAXPY( f, 1.0, f_utn );
	VecAXPY( f, 0.5*alpha, f_utnln );
	VecAXPY( f, 0.5*alpha, f_utnlp );
	VecAXPY( f, 0.5*alpha, f_utlp );
	VecAXPY( f, 1.0, pc->b_ut );
	pc->SolveLinSys( "tlc_", pc->velTop, pc->A1P, f, x );
	sol->UpdateField();
	VecDestroy( f );
	VecDestroy( x );
	delete sol;

	/* mass eqn corrector step */
	CreateVector( pc->etaInt, &f );
	CreateVector( pc->etaInt, &x );
	sol = new Vector( "eta-int-corr-sol", pc->etaInt, x, NULL, 0 );
	pc->AssembleMassEqnRHS( pc->velBotP, pc->etaIntP, f_dhup );
	VecAXPY( f, 1.0, f_ein );
	VecAXPY( f, 0.5, f_dhun );
	VecAXPY( f, 0.5, f_dhup );
	VecAXPY( f, 1.0, pc->b_ei );
	pc->SolveLinSys( "ihc_", pc->etaInt, pc->Ae, f, x );
	sol->UpdateField();
	VecDestroy( f );
	VecDestroy( x );
	delete sol;

	/* free the vectors */
	VecDestroy( f_utn );
	VecDestroy( f_ein );
	VecDestroy( f_dhun );
	VecDestroy( f_dhup );
	VecDestroy( f_utnln );
	VecDestroy( f_utnlp );
	VecDestroy( f_utlp );

	pc->velTopPrev->Copy( velTopTemp );
	delete velTopTemp;
}

void SolvePCFirstOrder( PredictorCorrector* pc ) {
	Vec 	f, x, f_utn, f_ein, f_dhun, f_utnln;
	Vector*	sol;

	/* vector steup */
	CreateVector( pc->velTop, &f_utn );		/* previous top layer velocity */
	CreateVector( pc->etaInt, &f_ein );		/* previous height interface perturbation */
	CreateVector( pc->etaInt, &f_dhun );	/* mass eqn nonlinear predictor */
	CreateVector( pc->velTop, &f_utnln );	/* top layer non linear predictor */

	pc->advTopN->Advect( pc->dt );

	/* top layer momentum eqn predictor step */
	CreateVector( pc->velTop, &f );
	CreateVector( pc->velTop, &x );
	sol = new Vector( "vel-top-pred-sol", pc->velTop, x, NULL, 0 );
	pc->AssemblePrevTimeStep( pc->advTopN->fieldSL, f_utn );
	pc->AssembleVelNonLinearRHS( pc->velTop, pc->etaInt, pc->tp, f_utnln );
	VecAXPY( f, 1.0, f_utn );
	VecAXPY( f, 1.0, f_utnln );
	VecAXPY( f, 1.0, pc->b_ut );
	pc->SolveLinSys( "tlp_", pc->velTop, pc->A1N, f, x );
	sol->UpdateField();
	VecDestroy( f );
	VecDestroy( x );
	delete sol;

	/* mass eqn predictor step */
	CreateVector( pc->etaInt, &f );
	CreateVector( pc->etaInt, &x );
	sol = new Vector( "eta-int-pred-sol", pc->etaInt, x, NULL, 0 );
	pc->AssemblePrevTimeStep( pc->etaInt, f_ein );
	pc->AssembleMassEqnRHS( pc->velBot, pc->etaInt, f_dhun );
	VecAXPY( f, 1.0, f_ein );
	VecAXPY( f, 1.0, f_dhun );
	VecAXPY( f, 1.0, pc->b_ei );
	pc->SolveLinSys( "ihp_", pc->etaInt, pc->Ae, f, x );
	sol->UpdateField();
	VecDestroy( f );
	VecDestroy( x );
	delete sol;

	/* free the vectors */
	VecDestroy( f_utn );
	VecDestroy( f_ein );
	VecDestroy( f_dhun );
	VecDestroy( f_utnln );

	pc->velTopPrev->Copy( pc->velTop );
}

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	bool			periodic[2]	= { true, false };
	bool			horiz[2]	= { false, false };
	bool			vert[2]		= { false, true };
	int			nx[2]		= { 16, 4 };
	double			min[2]		= { -8.0, +0.0 };
	double			max[2]		= { +8.0, +4.0 };
	Mesh*			vMesh		= new Mesh( "v-mesh", nx, "legendre", 8, min, max, periodic );
	Mesh*			hMesh		= new Mesh( "h-mesh", nx, "legendre", 6, min, max, periodic );
	BCs*			vbcs		= new BCs( vert, vert, horiz, horiz, vMesh, 0 );
	BCs*			hbcs		= new BCs( false, true, false, false, hMesh, 0 );
	Field*			vel		= new Field( "vel", vMesh, 2, vbcs );
	Field*			velPrev		= new Field( "vel-prev", vMesh, 2, vbcs );
	Field*			eta		= new Field( "eta", hMesh, 1, hbcs );
	Field*			topo		= new Field( "topo", hMesh, 1, NULL );
	double			time		= 0.0;
	double			dt		= 0.0075;
	SWParams		params;
	PredictorCorrector*	pc;

	params.L = params.H = params.H_i = params.Htilde = params.U = params.g = params.gPrime = params.rho = params.rho_0 = params.beta = 1.0;
	params.gamma = params.tau = params.f0 = params.nu = 0.0;

	PetscInitialize( &argc, &argv, (char)0, tag );

	vel->SetICFunc( 0, u0 );
	vel->SetICFunc( 1, v0 );
	eta->SetICFunc( 0, h0 );

	vMesh->Save();
	hMesh->Save();

	WriteXDMFHeader( 0 );
	WriteXDMF( &vel, 1, 0, time, dt );
	WriteXDMF( &eta, 1, 0, time, dt );
	WriteXDMFFooter( 0 );

	pc = new PredictorCorrector( vel, vel, eta, topo, &params, &params, dt, velPrev, velPrev );
	pc->unull = true;

	time += dt;
	pc->Assemble( 1 );
	SolvePCFirstOrder( pc );
	WriteXDMFHeader( 1 );
	WriteXDMF( &vel, 1, 1, time, dt );
	WriteXDMF( &eta, 1, 1, time, dt );
	WriteXDMFFooter( 1 );

	pc->Assemble( 2 );
	for( int step_i = 2; step_i <= 100; step_i++ ) {
		time += dt;
		SolvePC( pc );
		WriteXDMFHeader( step_i );
		WriteXDMF( &vel, 1, step_i, time, dt );
		WriteXDMF( &eta, 1, step_i, time, dt );
		WriteXDMFFooter( step_i );
	}

	delete vMesh;
	delete hMesh;
	delete vbcs;
	delete hbcs;
	delete vel;
	delete eta;
	delete topo;
	delete velPrev;
	delete pc;

	PetscFinalize();

	return EXIT_SUCCESS;
}
