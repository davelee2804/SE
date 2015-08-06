#include <string>
#include <iostream>

#include <hdf5.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SW.h"
#include "Params.h"
#include "2LRLOps.h"
#include "TwoLayerSW.h"

using namespace std;
using std::string;

TwoLayerSW::TwoLayerSW( Field* _velTop, Field* _velBot, Field* _etaInt, Field* _presSurf, TwoLayerParams* _params, double _dt ) {
	velTop		= _velTop;
	velBot		= _velBot;
	etaInt		= _etaInt;
	presSurf	= _presSurf;
	params		= _params;
	dt		= _dt;

	ass1 = false;
	ass2 = false;
	pnull = false;
	enull = false;
}

TwoLayerSW::~TwoLayerSW() {
	MatDestroy( PS );
	MatDestroy( PI );
	MatDestroy( UI );
	MatDestroy( VI );
}

void TwoLayerSW::Solve( Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev ) {
	Field*	velTopTemp	= NULL;
	Field*	velBotTemp	= NULL;
	Field*	etaIntTemp	= NULL;

	if( velTopPrev != NULL && velBotPrev != NULL && etaIntPrev != NULL ) {
		if( ass2 == false ) {
			MatDestroy( PS );
			MatDestroy( PI );
			MatDestroy( UI );
			MatDestroy( VI );

			AssembleMatrices( 1.5 );
			ass2 = true;
		}
		velTopTemp = new Field( "vel-top-temp", velTop->mesh, 2, NULL );
		velBotTemp = new Field( "vel-bot-temp", velBot->mesh, 2, NULL );
		etaIntTemp = new Field( "eta-int-temp", etaInt->mesh, 1, NULL );
		velTopTemp->Copy( velTop );
		velBotTemp->Copy( velBot );
		etaIntTemp->Copy( etaInt );

		SolveSecondOrder( velTopPrev, velBotPrev, etaIntPrev );

		velTopPrev->Copy( velTopTemp );
		velBotPrev->Copy( velBotTemp );
		etaIntPrev->Copy( etaIntTemp );
		delete velTopTemp;
		delete velBotTemp;
		delete etaIntTemp;
	}
	else {
		if( ass1 == false ) {
			AssembleMatrices( 1.0 );
			ass1 = true;
		}
		SolveFirstOrder();
	}
}

void TwoLayerSW::AssembleMatrices( double a ) {
	Laplacian*	lapMat		= new Laplacian( "laplacian", presSurf, presSurf, -dt*params->p/a );
	MassMatrix*	massMat		= new MassMatrix( "phi-mass-matrix", etaInt, etaInt, a );
	PhiMatrix*	phiMat		= new PhiMatrix( "phi-matrix", etaInt, etaInt, -0.5*dt*dt*params->g*params->H2, a, dt, params->f, params->beta );
	MassMatrix*	uMassMat	= new MassMatrix( "u-mass-matrix", velTop, velTop, a );
	BetaMatrix*	uBetaMat	= new BetaMatrix( "u-beta-matrix", velTop, velTop, dt, params->f, params->beta );
	StressTensor*	uViscMat;
	Operator**	ops;
	Vec		x;
	Vector*		sol;
	Matrix*		mat;

	/* surface pressure matrix */
	CreateVector( presSurf, &x );
	CreateMatrix( presSurf, presSurf, &PS );

	ops = new Operator*[1];
	ops[0] = lapMat;
	sol = new Vector( "p-sol", presSurf, x, NULL, 0 );
	mat = new Matrix( "p-mat", PS, sol, sol, NULL, ops, 1 );
	mat->Assemble();
	MatAssemblyBegin( PS, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( PS, MAT_FINAL_ASSEMBLY );

	delete sol;
	delete mat;
	VecDestroy( x );

	/* interfacial pressure matrix */
	CreateVector( etaInt, &x );
	CreateMatrix( etaInt, etaInt, &PI );

	sol = new Vector( "e-sol", etaInt, x, NULL, 0 );
	ops = new Operator*[2];
	ops[0] = massMat;
	ops[1] = phiMat;
	mat = new Matrix( "e-mat", PI, sol, sol, NULL, ops, 2 );
	mat->Assemble();
	MatAssemblyBegin( PI, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( PI, MAT_FINAL_ASSEMBLY );
	delete sol;
	delete mat;
	VecDestroy( x );

	/* create u matrix: TODO: add the boundary condition vector */
	CreateVector( velTop, &x );
	CreateMatrix( velTop, velTop, &UI );
	sol = new Vector( "u-sol", velTop, x, NULL, 0 );
	ops = new Operator*[2];
	ops[0] = uMassMat;
	ops[1] = uBetaMat;
	mat = new Matrix( "u-mat", UI, sol, sol, NULL, ops, 2 );
	mat->Assemble();
	MatAssemblyBegin( UI, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( UI, MAT_FINAL_ASSEMBLY );
	delete sol;
	delete mat;
	VecDestroy( x );

	/* create the baroclinic velocity matrix */
	uMassMat = new MassMatrix( "u-mass-matrix", velTop, velTop, a );
	uBetaMat = new BetaMatrix( "u-beta-matrix", velTop, velTop, dt, params->f, params->beta );
	uViscMat = new StressTensor( "u-visc-matrix", velTop, velTop, dt*params->nu );
	CreateVector( velTop, &x );
	CreateMatrix( velTop, velTop, &VI );
	sol = new Vector( "v-sol", velTop, x, NULL, 0 );
	ops = new Operator*[3];
	ops[0] = uMassMat;
	ops[1] = uBetaMat;
	ops[2] = uViscMat;
	mat = new Matrix( "v-mat", VI, sol, sol, NULL, ops, 3 );
	mat->Assemble();
	MatAssemblyBegin( VI, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( VI, MAT_FINAL_ASSEMBLY );
	delete sol;
	delete mat;
	VecDestroy( x );
}

void TwoLayerSW::CalcVelBar( Field* vel1, Field* vel2, Field* eta, Field* velb ) {
	double	*v1, *v2, e, h1, h2, htInv;

	htInv = 1.0/(params->H1 + params->H2);

	for( int i = 0; i < velb->mesh->nVertsTotal; i++ ) {
		v1 = vel1->vals[i];
		v2 = vel2->vals[i];
		eta->InterpGlobal( velb->mesh->verts[i], &e );
		h1 = params->H1 - e;
		h2 = params->H1 + e;

		velb->vals[i][0] = (h1*v1[0] + h2*v2[0])*htInv;
		velb->vals[i][1] = (h1*v1[1] + h2*v2[1])*htInv;
	}
}

void TwoLayerSW::AssembleGFirstOrder( Vec* G ) {
	Vector*		rhs;
	RHSOp**		ops		= new RHSOp*[1];
	GLayerVector2*	gCurrRHS	= NULL;

	CreateVector( velTop, G );

	gCurrRHS = new GLayerVector2( "g-curr", velTop->mesh, -dt, etaInt, velTop, velBot, 
				      params->H1, params->H2, params->nu, params->g, params->tau, params->k );
	ops[0] = gCurrRHS;
	rhs = new Vector( "G", velTop, *G, ops, 1 );
	rhs->Assemble();

	delete rhs;
}

void TwoLayerSW::AssembleGSecondOrder( Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev, Vec* G ) {
	Vector*		rhs;
	RHSOp**		ops		= new RHSOp*[2];
	GLayerVector2*	gCurrRHS;
	GLayerVector2*	gPrevRHS;

	CreateVector( velTop, G );

	gCurrRHS = new GLayerVector2( "g-curr", velTop->mesh, -2.0*dt, etaInt, velTop, velBot, 
				      params->H1, params->H2, params->nu, params->g, params->tau, params->k );
	gPrevRHS = new GLayerVector2( "g-prev", velTop->mesh, dt, etaIntPrev, velTopPrev, velBotPrev, 
				      params->H1, params->H2, params->nu, params->g, params->tau, params->k );
	ops[0] = gCurrRHS;
	ops[1] = gPrevRHS;
	rhs = new Vector( "G", velTop, *G, ops, 2 );
	rhs->Assemble();

	delete rhs;
}

void TwoLayerSW::SolvePressureFirstOrder() {
	Field*		velBar		= new Field( "vel-bar", velTop->mesh, 2, velTop->bcs );
	Field*		velHat		= new Field( "vel-hat", velTop->mesh, 2, velTop->bcs );
	Field*		deltaPres	= new Field( "delta-pres", presSurf->mesh, 1, presSurf->bcs );
	FieldRHS*	velRHS;
	GradFieldRHS*	gPresRHS;
	DivVelRHS*	divVelRHS;
	RHSOp**		ops;
	Vector*		rhs;
	Vector*		sol;
	Vec		x;
	Vec		f;
	Vec		G;
	KSP		ksp;
	MatNullSpace	null;

	CalcVelBar( velTop, velBot, etaInt, velBar );

	CreateVector( velBar, &f );
	CreateVector( velBar, &x );

	velRHS = new FieldRHS( "vel-rhs", velHat->mesh, 1.0, velBar );
	gPresRHS = new GradFieldRHS( "grad-pres-rhs", velHat->mesh, -dt*params->p, presSurf );
	ops = new RHSOp*[2];
	ops[0] = velRHS;
	ops[1] = gPresRHS;
	rhs = new Vector( "rhs", velHat, f, ops, 2 );
	sol = new Vector( "sol", velHat, x, NULL, 0 );
	rhs->Assemble();

	AssembleGFirstOrder( &G );
	VecAXPY( f, 1.0, G );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, UI, UI, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "vel-bar_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	KSPDestroy( ksp );

	VecDestroy( f );
	VecDestroy( x );
	VecDestroy( G );
	delete rhs;
	delete sol;
	delete velBar;

	CreateVector( presSurf, &f );
	CreateVector( presSurf, &x );

	divVelRHS = new DivVelRHS( "div-vel", presSurf->mesh, 1.0, velHat );
	ops = new RHSOp*[1];
	ops[0] = divVelRHS;
	rhs = new Vector( "rhs", deltaPres, f, ops, 1 );
	rhs->Assemble();

	sol = new Vector( "sol", deltaPres, x, NULL, 0 );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, PS, PS, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "pres-surf_" );
	KSPSetFromOptions( ksp );
	if( pnull ) {
		MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
		KSPSetNullSpace( ksp, null );
	}
	KSPSolve( ksp, f, x );
	if( pnull ) {
		MatNullSpaceDestroy( null );
	}
	KSPDestroy( ksp );
	sol->UpdateField();

	for( int i = 0; i < presSurf->mesh->nVertsTotal; i++ ) {
		presSurf->vals[i][0] += deltaPres->vals[i][0];
	}

	VecDestroy( f );
	VecDestroy( x );
	delete rhs;
	delete sol;
	delete velHat;
	delete deltaPres;
}

void TwoLayerSW::SolvePressureSecondOrder( Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev ) {
	Field*		velBarCurr	= new Field( "vel-bar-curr", velTop->mesh, 2, velTop->bcs );
	Field*		velBarPrev	= new Field( "vel-bar-prev", velTop->mesh, 2, velTop->bcs );
	Field*		velHat		= new Field( "vel-hat", velTop->mesh, 2, velTop->bcs );
	Field*		deltaPres	= new Field( "delta-pres", presSurf->mesh, 1, presSurf->bcs );
	FieldRHS*	velCurrRHS;
	FieldRHS*	velPrevRHS;
	GradFieldRHS*	gPresRHS;
	DivVelRHS*	divVelRHS;
	RHSOp**		ops;
	Vector*		rhs;
	Vector*		sol;
	Vec		x;
	Vec		f;
	Vec		G;
	KSP		ksp;
	MatNullSpace	null;

	CalcVelBar( velTop, velBot, etaInt, velBarCurr );
	CalcVelBar( velTopPrev, velBotPrev, etaIntPrev, velBarPrev );

	CreateVector( velBarCurr, &f );
	CreateVector( velBarCurr, &x );

	velCurrRHS = new FieldRHS( "vel-rhs", velHat->mesh, +2.0, velBarCurr );
	velPrevRHS = new FieldRHS( "vel-rhs", velHat->mesh, -0.5, velBarPrev );
	gPresRHS = new GradFieldRHS( "grad-pres-rhs", velHat->mesh, -dt*params->p, presSurf );
	ops = new RHSOp*[3];
	ops[0] = velCurrRHS;
	ops[1] = velPrevRHS;
	ops[2] = gPresRHS;
	rhs = new Vector( "rhs", velHat, f, ops, 3 );
	sol = new Vector( "sol", velHat, x, NULL, 0 );
	rhs->Assemble();

	AssembleGSecondOrder( velTopPrev, velBotPrev, etaIntPrev, &G );
	VecAXPY( f, 1.0, G );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, UI, UI, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "vel-bar_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	KSPDestroy( ksp );

	VecDestroy( f );
	VecDestroy( x );
	VecDestroy( G );
	delete rhs;
	delete sol;
	delete velBarCurr;
	delete velBarPrev;

	CreateVector( presSurf, &f );
	CreateVector( presSurf, &x );

	divVelRHS = new DivVelRHS( "div-vel", presSurf->mesh, 1.0, velHat );
	ops = new RHSOp*[1];
	ops[0] = divVelRHS;
	rhs = new Vector( "rhs", deltaPres, f, ops, 1 );
	rhs->Assemble();

	sol = new Vector( "sol", deltaPres, x, NULL, 0 );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, PS, PS, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "pres-surf_" );
	KSPSetFromOptions( ksp );
	if( pnull ) {
		MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
		KSPSetNullSpace( ksp, null );
	}
	KSPSolve( ksp, f, x );
	if( pnull ) {
		MatNullSpaceDestroy( null );
	}
	KSPDestroy( ksp );
	sol->UpdateField();

	for( int i = 0; i < presSurf->mesh->nVertsTotal; i++ ) {
		presSurf->vals[i][0] += deltaPres->vals[i][0];
	}

	VecDestroy( f );
	VecDestroy( x );
	delete rhs;
	delete sol;
	delete velHat;
	delete deltaPres;
}

void TwoLayerSW::SolveFirstOrder() {
	Advector*		velTopAdv	= new Advector( velTop, velTop );
	Advector*		velBotAdv	= new Advector( velBot, velBot );
	/* interfacial pressure rhs */
	FieldRHS*		etaRHS;
	DivHeightVelRHS*	dEtaVel1RHS;
	DivHeightVelRHS*	dEtaVel2RHS;
	PhiRHS*			u1RHS;
	PhiRHS*			u2RHS;
	/* momentum eqn rhs */
	FieldRHS*		velRHS;
	GradFieldRHS*		gPresRHS;
	WindStress2RHS*		windRHS;
	GradFieldRHS*		gEtaRHS;
	RHSOp**			ops;
	Vector*			sol;
	Vector*			rhs;
	Vec			f;
	Vec			x;
	KSP			ksp;
	MatNullSpace		null;

	//SolvePressureFirstOrder();

	CreateVector( etaInt, &f );
	CreateVector( etaInt, &x );

	velTopAdv->Advect( dt );
	velBotAdv->Advect( dt );

	sol = new Vector( "e-sol", etaInt, x, NULL, 0 );

	etaRHS      = new FieldRHS( "eta-rhs", etaInt->mesh, 1.0, etaInt );
	dEtaVel1RHS = new DivHeightVelRHS( "div-eta-vel-1-rhs", etaInt->mesh, -0.5*dt, etaInt, velTop );
	dEtaVel2RHS = new DivHeightVelRHS( "div-eta-vel-2-rhs", etaInt->mesh, -0.5*dt, etaInt, velBot );
	u1RHS       = new PhiRHS( "u1-rhs", etaInt->mesh, 1.0, velTop, +1.0, dt, params->H1, params->tau, params->k, 
				  params->nu, params->p, params->f, params->beta, presSurf, etaInt, velTopAdv->fieldSL );
	u2RHS       = new PhiRHS( "u2-rhs", etaInt->mesh, 1.0, velBot, -1.0, dt, params->H2, 0.0, 0.0, 
				  params->nu, params->p, params->f, params->beta, presSurf, etaInt, velBotAdv->fieldSL );
	ops = new RHSOp*[5];
	ops[0] = etaRHS;
	ops[1] = dEtaVel1RHS;
	ops[2] = dEtaVel2RHS;
	ops[3] = u1RHS;
	ops[4] = u2RHS;
	rhs = new Vector( "e-rhs", etaInt, f, ops, 5 );
	rhs->Assemble();
	
	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, PI, PI, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "eta-int_" );
	KSPSetFromOptions( ksp );
	if( enull ) {
		MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
		KSPSetNullSpace( ksp, null );
	}
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	if( enull ) {
		MatNullSpaceDestroy( null );
	}
	KSPDestroy( ksp );

	delete sol;
	delete rhs;
	VecDestroy( f );
	VecDestroy( x );

	/* update the upper layer velocity field */
	CreateVector( velTop, &x );
	CreateVector( velTop, &f );

	velRHS   = new FieldRHS( "vel-1-rhs", velTop->mesh, 1.0, velTopAdv->fieldSL );
	gPresRHS = new GradFieldRHS( "g-pres-rhs", velTop->mesh, -dt*params->p, presSurf );
	windRHS  = new WindStress2RHS( "wind-rhs", velTop->mesh, dt*params->tau, etaInt, params->k, params->H1 );
	ops = new RHSOp*[3];
	ops[0] = velRHS;
	ops[1] = gPresRHS;
	ops[2] = windRHS;
	rhs = new Vector( "u1-rhs", velTop, f, ops, 3 );
	sol = new Vector( "u1-sol", velTop, x, NULL, 0 );
	rhs->Assemble();

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, VI, VI, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "vel-top_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	KSPDestroy( ksp );

	delete sol;
	delete rhs;
	VecDestroy( f );
	VecDestroy( x );
	
	/* update the lower layer velocity field */
	CreateVector( velBot, &x );
	CreateVector( velBot, &f );

	velRHS   = new FieldRHS( "vel-2-rhs", velBot->mesh, 1.0, velBotAdv->fieldSL );
	gPresRHS = new GradFieldRHS( "g-pres-rhs", velBot->mesh, -dt*params->p, presSurf );
	gEtaRHS  = new GradFieldRHS( "g-eta-rhs", velBot->mesh, -dt*params->g, etaInt );
	ops = new RHSOp*[3];
	ops[0] = velRHS;
	ops[1] = gPresRHS;
	ops[2] = gEtaRHS;
	rhs = new Vector( "u2-rhs", velBot, f, ops, 3 );
	sol = new Vector( "u2-sol", velBot, x, NULL, 0 );
	rhs->Assemble();

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, VI, VI, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "vel-bot_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	KSPDestroy( ksp );

	delete sol;
	delete rhs;
	VecDestroy( f );
	VecDestroy( x );

	delete velTopAdv;
	delete velBotAdv;
}

void TwoLayerSW::SolveSecondOrder( Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev ) {
	Advector*		velTopAdv	= new Advector( velTop, velTop, velTopPrev, velTopPrev );
	Advector*		velBotAdv	= new Advector( velBot, velBot, velBotPrev, velBotPrev );
	/* interfacial pressure rhs ops */
	FieldRHS*		eta1RHS;
	FieldRHS*		eta2RHS;
	DivHeightVelRHS*	dEtaVel1CurrRHS;
	DivHeightVelRHS*	dEtaVel2CurrRHS;
	DivHeightVelRHS*	dEtaVel1PrevRHS;
	DivHeightVelRHS*	dEtaVel2PrevRHS;
	PhiRHS*			u1RHS;
	PhiRHS*			u2RHS;
	/* momentum eqn rhs ops */
	FieldRHS*		velRHS;
	GradFieldRHS*		gPresRHS;
	WindStress2RHS*		windRHS;
	GradFieldRHS*		gEtaRHS;
	RHSOp**			ops;
	Vector*			sol;
	Vector*			rhs;
	Vec			f;
	Vec			x;
	KSP			ksp;
	MatNullSpace		null;
	Field*			etaIntTemp	= new Field( "etaIntTemp", etaInt->mesh, 1, etaInt->bcs );
	Field*			velTopNew	= new Field( "vel-top-new", velTop->mesh, 2, velTop->bcs );
	Field*			velBotNew	= new Field( "vel-bot-new", velBot->mesh, 2, velBot->bcs );

	etaIntTemp->Copy( etaInt );
	for( int i = 0; i < velTop->mesh->nVertsTotal; i++ ) {
		velTopNew->vals[i][0] = 2.0*velTop->vals[i][0] - velTopPrev->vals[i][0];
		velTopNew->vals[i][1] = 2.0*velTop->vals[i][1] - velTopPrev->vals[i][1];
		velBotNew->vals[i][0] = 2.0*velBot->vals[i][0] - velBotPrev->vals[i][0];
		velBotNew->vals[i][1] = 2.0*velBot->vals[i][1] - velBotPrev->vals[i][1];
	}

	//SolvePressureSecondOrder( velTopPrev, velBotPrev, etaIntPrev );

	CreateVector( etaInt, &f );
	CreateVector( etaInt, &x );

	velTopAdv->Advect( dt );
	velBotAdv->Advect( dt );

	sol = new Vector( "e-sol", etaInt, x, NULL, 0 );

	eta1RHS         = new FieldRHS( "eta-curr-rhs", etaInt->mesh, +2.0, etaInt );
	eta2RHS         = new FieldRHS( "eta-prev-rhs", etaInt->mesh, -0.5, etaIntPrev );
	dEtaVel1CurrRHS = new DivHeightVelRHS( "d-eta-vel-1-curr-rhs", etaInt->mesh, -1.0*dt, etaInt, velTop );
	dEtaVel1PrevRHS = new DivHeightVelRHS( "d-eta-vel-1-prev-rhs", etaInt->mesh, +0.5*dt, etaIntPrev, velTopPrev );
	dEtaVel2CurrRHS = new DivHeightVelRHS( "d-eta-vel-2-curr-rhs", etaInt->mesh, -1.0*dt, etaInt, velBot );
	dEtaVel2PrevRHS = new DivHeightVelRHS( "d-eta-vel-2-prev-rhs", etaInt->mesh, +0.5*dt, etaIntPrev, velBotPrev );
	u1RHS           = new PhiRHS( "u1-rhs", etaInt->mesh, 1.5, velTopNew, +1.0, dt, params->H1, params->tau, params->k, 
	 			      params->nu, params->p, params->f, params->beta, presSurf, etaInt, velTopAdv->fieldSL );
	u2RHS           = new PhiRHS( "u2-rhs", etaInt->mesh, 1.5, velBotNew, -1.0, dt, params->H2, 0.0, 0.0, 
				      params->nu, params->p, params->f, params->beta, presSurf, etaInt, velBotAdv->fieldSL );
	ops = new RHSOp*[8];
	ops[0] = eta1RHS;
	ops[1] = eta2RHS;
	ops[2] = dEtaVel1CurrRHS;
	ops[3] = dEtaVel1PrevRHS;
	ops[4] = dEtaVel2CurrRHS;
	ops[5] = dEtaVel2PrevRHS;
	ops[6] = u1RHS;
	ops[7] = u2RHS;
	rhs = new Vector( "e-rhs", etaInt, f, ops, 8 );
	rhs->Assemble();
	
	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, PI, PI, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "eta-int_" );
	KSPSetFromOptions( ksp );
	if( enull ) {
		MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
		KSPSetNullSpace( ksp, null );
	}
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	if( enull ) {
		MatNullSpaceDestroy( null );
	}
	KSPDestroy( ksp );

	delete sol;
	delete rhs;
	VecDestroy( f );
	VecDestroy( x );

	/* update the upper layer velocity field */
	CreateVector( velTop, &x );
	CreateVector( velTop, &f );

	velRHS   = new FieldRHS( "vel-1-rhs", velTop->mesh, 1.0, velTopAdv->fieldSL );
	gPresRHS = new GradFieldRHS( "g-pres-rhs", velTop->mesh, -dt*params->p, presSurf );
	windRHS  = new WindStress2RHS( "wind-rhs", velTop->mesh, dt*params->tau, etaIntTemp, params->k, params->H1 );
	ops = new RHSOp*[3];
	ops[0] = velRHS;
	ops[1] = gPresRHS;
	ops[2] = windRHS;
	rhs = new Vector( "u1-rhs", velTop, f, ops, 3 );
	sol = new Vector( "u1-sol", velTop, x, NULL, 0 );
	rhs->Assemble();

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, VI, VI, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "vel-top_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	KSPDestroy( ksp );

	delete sol;
	delete rhs;
	VecDestroy( f );
	VecDestroy( x );

	/* update the lower layer velocity field */
	CreateVector( velBot, &x );
	CreateVector( velBot, &f );

	velRHS   = new FieldRHS( "vel-2-rhs", velBot->mesh, 1.0, velBotAdv->fieldSL );
	gPresRHS = new GradFieldRHS( "g-pres-rhs", velBot->mesh, -dt*params->p, presSurf );
	gEtaRHS  = new GradFieldRHS( "g-eta-rhs", velBot->mesh, -dt*params->g, etaInt );
	ops = new RHSOp*[3];
	ops[0] = velRHS;
	ops[1] = gPresRHS;
	ops[2] = gEtaRHS;
	rhs = new Vector( "u2-rhs", velBot, f, ops, 3 );
	sol = new Vector( "u2-sol", velBot, x, NULL, 0 );
	rhs->Assemble();

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, VI, VI, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "vel-bot_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	sol->UpdateField();
	KSPDestroy( ksp );

	delete sol;
	delete rhs;
	VecDestroy( f );
	VecDestroy( x );

	delete velTopAdv;
	delete velBotAdv;
	delete etaIntTemp;
	delete velTopNew;
	delete velBotNew;
}
