#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SW.h"
#include "Params.h"
#include "2LRLOps.h"
#include "BarotropicEqn.h"

using namespace std;
using std::string;

BarotropicEqn::BarotropicEqn( Field* _velocity, Field* _pressure, double _dt, Params* _params, Field* _topo ) {
	double t = 0;

	velocity 	= _velocity;
	pressure 	= _pressure;
	dt       	= _dt;
	params 		= _params;
	topo		= _topo;

	totalHeight = new Field( "total-height", velocity->mesh, 1, NULL );
	for( int node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		if( topo != NULL ) { 
			topo->InterpGlobal( velocity->mesh->verts[node_i], &t ); 
		}
		totalHeight->vals[node_i][0] = params->H - t;
	}

	M_u = PETSC_NULL;
	M_p = PETSC_NULL;
	b_u = PETSC_NULL;

	pnull = false;
}

BarotropicEqn::~BarotropicEqn() {
	if( M_u != PETSC_NULL ) MatDestroy( M_u );
	if( M_p != PETSC_NULL ) MatDestroy( M_p );
	if( b_u != PETSC_NULL ) VecDestroy( b_u );

	delete totalHeight;
}

//#define TEST

void BarotropicEqn::Setup( int order ) {
	MassMatrix*		massMat;
#ifndef TEST
	BetaMatrix*		betaMat;
#else
	StressTensor*		viscMat;
#endif
	//NonLinearLaplacian*	lapMat;
	Laplacian*		lapMat;
	SVV*			svvMat;
	Operator**		ops;
	Vector*			uSol;
	Vector*			pSol;
	Vector*			uRhs;
	Matrix*			uMat;
	Matrix*			pMat;
	Vec			u;
	Vec			p;
	double			alpha		= ( order == 2 ) ? 1.5 : 1.0;

	if( M_u != PETSC_NULL ) { MatDestroy( M_u ); }
	if( M_p != PETSC_NULL ) { MatDestroy( M_p ); }
	if( b_u != PETSC_NULL ) { VecDestroy( b_u ); }

	CreateMatrix( velocity, velocity, &M_u );
	CreateMatrix( pressure, pressure, &M_p );

	CreateVector( velocity, &u );
	CreateVector( pressure, &p );
	CreateVector( velocity, &b_u );

	uSol = new Vector( "u-sol", velocity, u, NULL, 0 );
	pSol = new Vector( "p-sol", pressure, p, NULL, 0 );
	uRhs = new Vector( "u-rhs", velocity, b_u, NULL, 0 );

	/* assemble the velocity step matrix */
	massMat = new MassMatrix( "u-mass-mat", velocity, velocity, alpha );
#ifndef TEST
	betaMat = new BetaMatrix( "u-beta-mat", velocity, velocity, dt, params->f0, params->beta );
#else
	viscMat = new StressTensor( "u-visc-mat", velocity, velocity, dt*params->nu );
#endif
	svvMat  = new SVV( "u-svv-mat", velocity, velocity, dt/velocity->mesh->el->N, params->svv );
	ops = new Operator*[3];
	ops[0] = massMat;
#ifndef TEST
	ops[1] = betaMat;
#else
	ops[1] = viscMat;
#endif
	ops[2] = svvMat;
	uMat = new Matrix( "M_u", M_u, uSol, uSol, uRhs, ops, 3 );
	uMat->Assemble();
	MatAssemblyBegin( M_u, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( M_u, MAT_FINAL_ASSEMBLY );

	/* assemble the pressure step matrix */
	//lapMat = new NonLinearLaplacian( "p-lap-mat", pressure, pressure, -dt*params->pFac, totalHeight );
	lapMat = new Laplacian( "p-lap-mat", pressure, pressure, -dt*params->pFac*params->H/alpha );
	ops = new Operator*[1];
	ops[0] = lapMat;
	pMat = new Matrix( "M_p", M_p, pSol, pSol, NULL, ops, 1 );
	pMat->Assemble();
	MatAssemblyBegin( M_p, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( M_p, MAT_FINAL_ASSEMBLY );

	delete uSol;
	delete pSol;
	delete uRhs;
	delete uMat;
	delete pMat;
	VecDestroy( u );
	VecDestroy( p );
}

void BarotropicEqn::Solve( Field* velPrev, Vec G ) {
	Field*			velHat		= new Field( "vel-hat", velocity->mesh, 2, velocity->bcs );
	Field*			presDelta	= new Field( "pres-delta", pressure->mesh, 1, pressure->bcs );
	FieldRHS*		velCurrRHS;
	FieldRHS*		velPrevRHS	= NULL;
	GradFieldRHS*		gradPresRHS;
	DivHeightVelRHS*	divHVelRHS;
	RHSOp**			ops;
	Vector*			uSol;
	Vector*			pSol;
	Vector*			uRhs;
	Vector*			pRhs;
	Vec			u;
	Vec			p;
	Vec			fu;
	Vec			fp;
	double			alpha		= ( velPrev ) ? 2.0/3.0 : 1.0;
	double			delta		= ( velPrev ) ? 2.0 : 1.0;
	int			nVelOps		= ( velPrev ) ? 3 : 2;
	KSP			ksp;
	double			**dp;
	MatNullSpace		null;
	//Vec			p_neumann;

	CreateVector( velocity, &u );
	CreateVector( pressure, &p );
	CreateVector( velocity, &fu );
	CreateVector( pressure, &fp );

	uSol = new Vector( "u-sol", velHat, u, NULL, 0 );
	pSol = new Vector( "p-sol", presDelta, p, NULL, 0 );

	/* assemble the velocity equation rhs vector */
	ops = new RHSOp*[nVelOps];
	velCurrRHS = new FieldRHS( "vel-curr-rhs", velocity->mesh, delta, velocity );
	gradPresRHS = new GradFieldRHS( "grad-pres-rhs", velocity->mesh, -dt*params->pFac, pressure );
	if( velPrev ) { velPrevRHS = new FieldRHS( "vel-prev-rhs", velocity->mesh, -0.5, velPrev ); }
	ops[0] = velCurrRHS;
	ops[1] = gradPresRHS;
	if( velPrev ) { ops[2] = velPrevRHS; }
	uRhs = new Vector( "v-rhs", velHat, fu, ops, nVelOps );
	uRhs->Assemble();

	VecAXPY( fu, 1.0, b_u ); /* ??? for periodic bcs? */
	VecAXPY( fu, 1.0, G );

	/* solve the velocity equation */
	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, M_u, M_u, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "bt_u_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, fu, u );
	KSPDestroy( ksp );
	uSol->UpdateField();

	/* assemble the pressure equation rhs vector */
	divHVelRHS = new DivHeightVelRHS( "div-height-vel-rhs", pressure->mesh, 1.0, totalHeight, velHat );
	ops = new RHSOp*[1];
	ops[0] = divHVelRHS;
	pRhs = new Vector( "p-rhs", presDelta, fp, ops, 1 );
	pRhs->Assemble();

	/* add the neumann bc */
	//AssemblePressureBC( velHat, &p_neumann );
	//VecAXPY( fp, 1.0, p_neumann );

	/* solve the pressure equation */
	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, M_p, M_p, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "bt_p_" );
	KSPSetFromOptions( ksp );
	if( pnull ) {
		MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
		KSPSetNullSpace( ksp, null );
	}
	KSPSolve( ksp, fp, p );
	KSPDestroy( ksp );
	if( pnull ) {
		MatNullSpaceDestroy( null );
	}
	pSol->UpdateField();

	/* delete the neumann bc */
	//VecDestroy( p_neumann );

	/* update the velocity */
	dp = new double*[1];
	dp[0] = new double[2];
	for( int node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		presDelta->InterpDerivsGlobal( velocity->mesh->verts[node_i], dp );
		if( !velocity->bcs->IsBCNodeAndDof( node_i, 0 ) ) {
			velocity->vals[node_i][0] = velHat->vals[node_i][0] - dt*params->pFac*alpha*dp[0][0];
		}
		if( !velocity->bcs->IsBCNodeAndDof( node_i, 1 ) ) {
			velocity->vals[node_i][1] = velHat->vals[node_i][1] - dt*params->pFac*alpha*dp[0][1];
		}
	}
	delete[] dp[0];
	delete[] dp;

	for( int node_i = 0; node_i < pressure->mesh->nVertsTotal; node_i++ ) {
		pressure->vals[node_i][0] += presDelta->vals[node_i][0];
	}

	delete velHat;
	delete presDelta;

	delete uSol;
	delete pSol;
	delete uRhs;
	delete pRhs;

	VecDestroy( u );
	VecDestroy( p );
	VecDestroy( fu );
	VecDestroy( fp );
}

void BarotropicEqn::AssembleGFirstOrder( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* etaInt, Vec* G ) {
	Vector*		rhs;
	RHSOp**		ops;
	Field*		hTop;
	Field*		hBot;
	GLayerVector*	gCurrRHS;

	CreateVector( velocity, G );
	CalcHeights( tp, bp, etaInt, &hTop, &hBot );

	gCurrRHS = new GLayerVector( "g-curr", velTop->mesh, -dt, topo, hTop, hBot, velTop, velBot, 
				     tp->H, bp->H, params->nu, params->gamma, params->hFac, params->tau, params->k );
	ops = new RHSOp*[1];
	ops[0] = gCurrRHS;
	rhs = new Vector( "G", velocity, *G, ops, 1 );
	rhs->Assemble();

	delete rhs;
	delete hTop;
	delete hBot;
}

void BarotropicEqn::AssembleGSecondOrder( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* etaInt, Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev, Vec* G ) {
	Vector*		rhs;
	RHSOp**		ops;
	Field*		hTop;
	Field*		hBot;
	Field*		hTopPrev;
	Field*		hBotPrev;
	GLayerVector*	gCurrRHS;
	GLayerVector*	gPrevRHS;

	CreateVector( velocity, G );

	CalcHeights( tp, bp, etaInt, &hTop, &hBot );
	CalcHeights( tp, bp, etaIntPrev, &hTopPrev, &hBotPrev );

	gCurrRHS = new GLayerVector( "g-curr", velTop->mesh, -2.0*dt, topo, hTop, hBot, velTop, velBot, 
				     tp->H, bp->H, params->nu, params->gamma, params->hFac, params->tau, params->k );
	gPrevRHS = new GLayerVector( "g-prev", velTop->mesh, dt, topo, hTopPrev, hBotPrev, velTopPrev, velBotPrev, 
				     tp->H, bp->H, params->nu, params->gamma, params->hFac, params->tau, params->k );
	ops = new RHSOp*[2];
	ops[0] = gCurrRHS;
	ops[1] = gPrevRHS;
	rhs = new Vector( "G", velocity, *G, ops, 2 );
	rhs->Assemble();

	delete rhs;
	delete hTop;
	delete hBot;
	delete hTopPrev;
	delete hBotPrev;
}

void BarotropicEqn::CalcVelBar( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* etaInt, Field* velBar ) {
	double	*v_1, *v_2, h_1, h_2, e, b = 0.0, hinv;

	for( int node_i = 0; node_i < velBar->mesh->nVertsTotal; node_i++ ) {
		v_1 = velTop->vals[node_i];
		v_2 = velBot->vals[node_i];
		etaInt->InterpGlobal( velBar->mesh->verts[node_i], &e );
		if( topo != NULL ) {
			topo->InterpGlobal( velBar->mesh->verts[node_i], &b );
		}
		h_1 = tp->H - e;
		h_2 = bp->H + e - b;
		hinv = 1.0/( h_1 + h_2 );
		velBar->vals[node_i][0] = hinv*( h_1*v_1[0] + h_2*v_2[0] );
		velBar->vals[node_i][1] = hinv*( h_1*v_1[1] + h_2*v_2[1] );
	}
}

void BarotropicEqn::CalcHeights( Params* tp, Params* bp, Field* eta, Field** height1, Field** height2 ) {
	double e, b = 0.0;

	*height1 = new Field( "height1", eta->mesh, 1, NULL );
	*height2 = new Field( "height2", eta->mesh, 1, NULL );

	for( int node_i = 0; node_i < eta->mesh->nVertsTotal; node_i++ ) {
		e = eta->vals[node_i][0];
		if( topo != NULL ) {
			topo->InterpGlobal( eta->mesh->verts[node_i], &b );
		}
		(*height1)->vals[node_i][0] = tp->H - e;
		(*height2)->vals[node_i][0] = bp->H + e - b;
	}
}
/*
void BarotropicEqn::AssemblePressureBC( Params* tp, Params* bp, Params* bt, 
					Field* velTop, Field* velTopPrev, 
					Field* velBot, Field* velBotPrev, 
					Field* etaInt, Field* etaIntPrev, 
					Field* velBar, Field* velBarPrev, Vec* n ) 
{
	Field			*height1, *height2, *height1Prev, *height2Prev;
	BaroclinicPressureBC*	pNeumannBCCurr;
	BaroclinicPressureBC*	pNeumannBCPrev;
	RHSOp**			ops		= new RHSOp*[2];
	Vector*			rhs;

	CreateVector( pressure, n );

	CalcHeights( tp, bp, etaInt, &height1, &height2 );
	CalcHeights( tp, bp, etaIntPrev, &height1Prev, &height2Prev );

	pNeumannBCCurr = new BaroclinicPressureBC( "neumann-pressure-rhs-curr", pressure->mesh, +2.0*dt/3.0, NULL, bt, velTop, velBot, velBar, height1, height2 );
	pNeumannBCPrev = new BaroclinicPressureBC( "heumann-pressure-rhs-prev", pressure->mesh, -2.0*dt/3.0, NULL, bt, velTopPrev, velBotPrev, velBarPrev, height1Prev, height2Prev );
	ops[0] = pNeumannBCCurr;
	ops[1] = pNeumannBCPrev;
	rhs = new Vector( "p-neumann-rhs", pressure, *n, ops, 2 );	
	rhs->Assemble();
	delete rhs;

	delete height1;
	delete height2;
	delete height1Prev;
	delete height2Prev;
}
*/
void BarotropicEqn::AssemblePressureBC( Field* velHat, Vec* n ) {
	BarotropicPressureBC*	pNeumannBC;
	RHSOp**			ops		= new RHSOp*[1];
	Vector*			rhs;

	CreateVector( pressure, n );

	pNeumannBC = new BarotropicPressureBC( "neumann-pressure-rhs", pressure->mesh, -1.0, velHat );
	ops[0] = pNeumannBC;
	rhs = new Vector( "p-neumann-rhs", pressure, *n, ops, 1 );	
	rhs->Assemble();
	delete rhs;
}
