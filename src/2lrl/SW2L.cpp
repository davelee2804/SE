#include <string>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SW.h"
#include "Params.h"
#include "2LRLOps.h"
#include "SW2L.h"

using namespace std;
using std::string;

SW2L::SW2L( Field* _velTop, Field* _velBot, Field* _etaInt, Field* _velTopPrev, Field* _velBotPrev, Field* _etaIntPrev, 
	    Field* _presSurf, Field* _bTopog, Params* _tp, Params* _bp, double _dt ) 
{
	int	size	= 0;

	velTop		= _velTop;
	velBot		= _velBot;
	etaInt		= _etaInt;
	velTopPrev	= _velTopPrev;
	velBotPrev	= _velBotPrev;
	etaIntPrev	= _etaIntPrev;
	presSurf	= _presSurf;
	bTopog		= _bTopog;
	tp		= _tp;
	bp		= _bp;
	dt		= _dt;

	velTopTemp = new Field( "vel-top-temp", velTop->mesh, 2, NULL );
	velBotTemp = new Field( "vel-bot-temp", velBot->mesh, 2, NULL );
	etaIntTemp = new Field( "eta-int-temp", etaInt->mesh, 1, NULL );

	size += ( 2*velTop->mesh->nVertsTotal - velTop->bcs->size[0] - velTop->bcs->size[1] );
	size += ( 2*velBot->mesh->nVertsTotal - velBot->bcs->size[0] - velBot->bcs->size[1] );
	size += (   etaInt->mesh->nVertsTotal - etaInt->bcs->size[0]                        );

	MatCreate( MPI_COMM_WORLD, &A );
	MatSetSizes( A, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( A, MATSEQAIJ );
	MatSeqAIJSetPreallocation( A, 2*4*2*velTop->mesh->el->nNodes + 4*etaInt->mesh->el->nNodes, PETSC_NULL );

	VecCreate( MPI_COMM_WORLD, &b );
	VecSetSizes( b, size, PETSC_DETERMINE );
	VecSetType( b, VECSEQ );
	VecSetOption( b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );

	unull = vnull = false;
	enull = true;
}

SW2L::~SW2L() {
	delete velTopTemp;
	delete velBotTemp;
	delete etaIntTemp;
	MatDestroy( A );
	VecDestroy( b );
}

void SW2L::Assemble( int order ) {
	MassMatrix*		massMat;
	BetaMatrix*		betaMat;
	StressTensor*		viscMat;
	SVV*			svvMat;
	Gradient*		gradMat;
	Divergence*		diveMat;
	Operator**		ops;
	Vec			x;
	Vector			*row, *col, *rhs;
	Matrix			*mat;
	double			a 		= ( order == 2 ) ? 1.5 : 1.0;

	MatZeroEntries( A );

	/* u1-u1 block */
	VecDuplicate( b, &x );
	row = new Vector( "row", velTop, x, NULL, 0 );
	col = new Vector( "col", velTop, x, NULL, 0 );
	rhs = new Vector( "rhs", velTop, b, NULL, 0 );
	massMat = new MassMatrix( "u1-u1-mass-mat", velTop, velTop, a + dt*tp->gamma );
	betaMat = new BetaMatrix( "u1-u1-beta-mat", velTop, velTop, dt, tp->f0, tp->beta );
	viscMat = new StressTensor( "u1-u1-visc-mat", velTop, velTop, dt*tp->nu );
	svvMat  = new SVV( "u1-u1-svv-mat", velTop, velTop, dt/velTop->mesh->el->N, tp->svv );
	ops = new Operator*[4];
	ops[0] = massMat;
	ops[1] = betaMat;
	ops[2] = viscMat;
	ops[3] = svvMat;
	mat = new Matrix( "A-u1-u1", A, row, col, rhs, ops, 4 );
	mat->Assemble();
	VecDestroy( x );
	delete row;
	delete col;
	delete rhs;
	delete mat;

	/* u2-u2 block */
	VecDuplicate( b, &x );
	row = new Vector( "row", velBot, x, NULL, 0 );
	col = new Vector( "col", velBot, x, NULL, 0 );
	rhs = new Vector( "rhs", velBot, b, NULL, 0 );
	massMat = new MassMatrix( "u2-u2-mass-mat", velBot, velBot, a + dt*bp->gamma );
	betaMat = new BetaMatrix( "u2-u2-beta-mat", velBot, velBot, dt, bp->f0, bp->beta );
	viscMat = new StressTensor( "u2-u2-visc-mat", velBot, velBot, dt*bp->nu );
	svvMat  = new SVV( "u2-u2-svv-mat", velBot, velBot, dt/velBot->mesh->el->N, bp->svv );
	ops = new Operator*[4];
	ops[0] = massMat;
	ops[1] = betaMat;
	ops[2] = viscMat;
	ops[3] = svvMat;
	mat = new Matrix( "A-u2-u2", A, row, col, rhs, ops, 4 );
	mat->Assemble();
	VecDestroy( x );
	delete row;
	delete col;
	delete rhs;
	delete mat;

	/* u2-ei block */
	VecDuplicate( b, &x );
	row = new Vector( "row", velBot, x, NULL, 0 );
	col = new Vector( "col", etaInt, x, NULL, 0 );
	rhs = new Vector( "rhs", velBot, b, NULL, 0 );
	gradMat = new Gradient( "u2-ei-grad-mat", velBot, etaInt, dt*bp->hFac );
	ops = new Operator*[1];
	ops[0] = gradMat;
	mat = new Matrix( "A-u2-ei", A, row, col, rhs, ops, 1 );
	mat->Assemble();
	VecDestroy( x );
	delete row;
	delete col;
	delete rhs;
	delete mat;

	/* ei-u1 block */
	VecDuplicate( b, &x );
	row = new Vector( "row", etaInt, x, NULL, 0 );
	col = new Vector( "col", velTop, x, NULL, 0 );
	rhs = new Vector( "rhs", etaInt, b, NULL, 0 );
	diveMat = new Divergence( "ei-u1-dive-mat", etaInt, velTop, -0.5*dt*tp->H );
	ops = new Operator*[1];
	ops[0] = diveMat;
	mat = new Matrix( "A-ei-u1", A, row, col, rhs, ops, 1 );
	mat->Assemble();
	VecDestroy( x );
	delete row;
	delete col;
	delete rhs;
	delete mat;

	/* ei-u2 block */
	VecDuplicate( b, &x );
	row = new Vector( "row", etaInt, x, NULL, 0 );
	col = new Vector( "col", velBot, x, NULL, 0 );
	rhs = new Vector( "rhs", etaInt, b, NULL, 0 );
	diveMat = new Divergence( "ei-u2-dive-mat", etaInt, velBot, +0.5*dt*bp->H );
	ops = new Operator*[1];
	ops[0] = diveMat;
	mat = new Matrix( "A-ei-u2", A, row, col, rhs, ops, 1 );
	mat->Assemble();
	VecDestroy( x );
	delete row;
	delete col;
	delete rhs;
	delete mat;

	/* ei-ei block */
	VecDuplicate( b, &x );
	row = new Vector( "row", etaInt, x, NULL, 0 );
	col = new Vector( "col", etaInt, x, NULL, 0 );
	rhs = new Vector( "rhs", etaInt, b, NULL, 0 );
	massMat = new MassMatrix( "ei-ei-mass-mat", etaInt, etaInt, a );
	ops = new Operator*[1];
	ops[0] = massMat;
	mat = new Matrix( "A-ei-ei", A, row, col, rhs, ops, 1 );
	mat->Assemble();
	VecDestroy( x );
	delete row;
	delete col;
	delete rhs;
	delete mat;

	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
}

void SW2L::Solve( int order ) {
	velTopTemp->Copy( velTop );
	velBotTemp->Copy( velBot );
	etaIntTemp->Copy( etaInt );

	if( order == 1 ) { FirstOrder(); }
	if( order == 2 ) { SecondOrder(); }

	velTopPrev->Copy( velTopTemp );
	velBotPrev->Copy( velBotTemp );
	etaIntPrev->Copy( etaIntTemp );
}

void SW2L::FirstOrder() {
	Advector*		topAdv		= new Advector( velTop, velTop );
	Advector*		botAdv		= new Advector( velBot, velBot );
	FieldRHS*		currRhs;
	GradFieldRHS*		gPresRhs;
	//GradFieldRHS*		gTopoRhs;
	WindStressRHS*		windRhs;
	DivHeightVelRHS*	dEtaVelTopRhs;
	DivHeightVelRHS*	dEtaVelBotRhs;
	RHSOp**			ops;
	Vec			x, f;
	Vector			*utSol, *ubSol, *eiSol, *rhs;

	VecDuplicate( b, &x );
	VecDuplicate( b, &f );
	VecZeroEntries( x );
	VecZeroEntries( f );

	topAdv->Advect( dt );
	botAdv->Advect( dt );

	utSol = new Vector( "ut-sol", velTop, x, NULL, 0 );
	ubSol = new Vector( "ub-sol", velBot, x, NULL, 0 );
	eiSol = new Vector( "ei-sol", etaInt, x, NULL, 0 );

	/* u1 rhs */
	currRhs    = new FieldRHS( "vel-top-adv-rhs", velTop->mesh, 1.0, topAdv->fieldSL );
	gPresRhs   = new GradFieldRHS( "grad-pres-rhs", velTop->mesh, -1.0*dt*tp->pFac, presSurf );
	windRhs    = new WindStressRHS( "wind-stress-rhs", velTop->mesh, dt*tp->tau, etaInt, tp->k, tp->H );
	ops = new RHSOp*[3];
	ops[0] = currRhs;
	ops[1] = gPresRhs;
	ops[2] = windRhs;
	rhs = new Vector( "u1-rhs", velTop, f, ops, 3 );
	rhs->Assemble();
	delete rhs;

	/* u2 rhs */
	currRhs    = new FieldRHS( "vel-bot-adv-rhs", velBot->mesh, 1.0, botAdv->fieldSL );
	gPresRhs   = new GradFieldRHS( "grad-pres-rhs", velBot->mesh, -1.0*dt*bp->pFac, presSurf );
	//gTopoRhs   = new GradFieldRHS( "grad-topo-rhs", velBot->mesh, -1.0*dt*bp->bFac, bTopog );
	ops = new RHSOp*[2];
	ops[0] = currRhs;
	ops[1] = gPresRhs;
	//ops[2] = gTopoRhs;
	rhs = new Vector( "u2-rhs", velBot, f, ops, 2 );
	rhs->Assemble();
	delete rhs;

	/* ei rhs */
	currRhs    = new FieldRHS( "eta-curr-rhs", etaInt->mesh, 1.0, etaInt );
	dEtaVelTopRhs = new DivHeightVelRHS( "div-eta-vel-top-rhs", etaInt->mesh, -0.5*dt, etaInt, velTop );
	dEtaVelBotRhs = new DivHeightVelRHS( "div-eta-vel-bot-rhs", etaInt->mesh, -0.5*dt, etaInt, velBot );
	ops = new RHSOp*[3];
	ops[0] = currRhs;
	ops[1] = dEtaVelTopRhs;
	ops[2] = dEtaVelBotRhs;
	rhs = new Vector( "ei-rhs", etaInt, f, ops, 3 );
	rhs->Assemble();
	delete rhs;

	VecAXPY( f, 1.0, b );
	SolveLinSys( f, x );
	utSol->UpdateField();
	ubSol->UpdateField();
	eiSol->UpdateField();

	delete topAdv;
	delete botAdv;
	delete utSol;
	delete ubSol;
	delete eiSol;
	VecDestroy( x );
	VecDestroy( f );
}

void SW2L::SecondOrder() {
	Advector*		topAdv		= new Advector( velTop, velTop, velTopPrev, velTopPrev );
	Advector*		botAdv		= new Advector( velBot, velBot, velBotPrev, velBotPrev );
	FieldRHS*		currRhs;
	FieldRHS*		prevRhs;
	GradFieldRHS*		gPresRhs;
	//GradFieldRHS*		gTopoRhs;
	WindStressRHS*		windRhs;
	WindStressRHS*		windPrevRhs;
	DivHeightVelRHS*	dEtaVelTopRhs;
	DivHeightVelRHS*	dEtaVelBotRhs;
	DivHeightVelRHS*	dEtaVelTopPrevRhs;
	DivHeightVelRHS*	dEtaVelBotPrevRhs;
	RHSOp**			ops;
	Vec			x, f;
	Vector			*utSol, *ubSol, *eiSol, *rhs;

	VecDuplicate( b, &x );
	VecDuplicate( b, &f );
	VecZeroEntries( x );
	VecZeroEntries( f );

	topAdv->Advect( dt );
	botAdv->Advect( dt );

	utSol = new Vector( "ut-sol", velTop, x, NULL, 0 );
	ubSol = new Vector( "ub-sol", velBot, x, NULL, 0 );
	eiSol = new Vector( "ei-sol", etaInt, x, NULL, 0 );

	/* u1 rhs */
	currRhs        = new FieldRHS( "vel-top-adv-rhs", velTop->mesh, 1.0, topAdv->fieldSL );
	gPresRhs       = new GradFieldRHS( "grad-pres-rhs", velTop->mesh, -1.0*dt*tp->pFac, presSurf );
	windRhs        = new WindStressRHS( "wind-stress-rhs", velTop->mesh, 2.0*dt*tp->tau, etaInt, tp->k, tp->H );
	windPrevRhs    = new WindStressRHS( "wind-stress-prev-rhs", velTop->mesh, -1.0*dt*tp->tau, etaIntPrev, tp->k, tp->H );
	ops = new RHSOp*[4];
	ops[0] = currRhs;
	ops[1] = gPresRhs;
	ops[2] = windRhs;
	ops[3] = windPrevRhs;
	rhs = new Vector( "u1-rhs", velTop, f, ops, 4 );
	rhs->Assemble();
	delete rhs;

	/* u2 rhs */
	currRhs        = new FieldRHS( "vel-bot-adv-rhs", velBot->mesh, 1.0, botAdv->fieldSL );
	gPresRhs       = new GradFieldRHS( "grad-pres-rhs", velBot->mesh, -1.0*dt*bp->pFac, presSurf );
	//gTopoRhs       = new GradFieldRHS( "grad-topo-rhs", velBot->mesh, -1.0*dt*bp->bFac, bTopog );
	ops = new RHSOp*[2];
	ops[0] = currRhs;
	ops[1] = gPresRhs;
	//ops[2] = gTopoRhs;
	rhs = new Vector( "u2-rhs", velBot, f, ops, 2 );
	rhs->Assemble();
	delete rhs;

	/* ei rhs */
	currRhs           = new FieldRHS( "eta-curr-rhs", etaInt->mesh, +2.0, etaInt );
	prevRhs           = new FieldRHS( "eta-prev-rhs", etaInt->mesh, -0.5, etaIntPrev );
	dEtaVelTopRhs     = new DivHeightVelRHS( "div-eta-vel-top-rhs", etaInt->mesh, -0.5*2.0*dt, etaInt, velTop );
	dEtaVelBotRhs     = new DivHeightVelRHS( "div-eta-vel-bot-rhs", etaInt->mesh, -0.5*2.0*dt, etaInt, velBot );
	dEtaVelTopPrevRhs = new DivHeightVelRHS( "div-eta-vel-top-prev-rhs", etaInt->mesh, +0.5*dt, etaIntPrev, velTopPrev );
	dEtaVelBotPrevRhs = new DivHeightVelRHS( "div-eta-vel-bot-prev-rhs", etaInt->mesh, +0.5*dt, etaIntPrev, velBotPrev );
	ops = new RHSOp*[6];
	ops[0] = currRhs;
	ops[1] = prevRhs;
	ops[2] = dEtaVelTopRhs;
	ops[3] = dEtaVelBotRhs;
	ops[4] = dEtaVelTopPrevRhs;
	ops[5] = dEtaVelBotPrevRhs;
	rhs = new Vector( "ei-rhs", etaInt, f, ops, 6 );
	rhs->Assemble();
	delete rhs;

	VecAXPY( f, 1.0, b );
	SolveLinSys( f, x );
	utSol->UpdateField();
	ubSol->UpdateField();
	eiSol->UpdateField();

	delete topAdv;
	delete botAdv;
	delete utSol;
	delete ubSol;
	delete eiSol;
	VecDestroy( x );
	VecDestroy( f );
}

void SW2L::SolveLinSys( Vec f, Vec x ) {
	KSP			ksp;
	KSP*			subksp;
	PC			pc;
	MatNullSpace		null;
	IS			utis, vtis, ubis, vbis, eis;
	int			uSize		= velTop->mesh->nVertsTotal - velTop->bcs->size[0];
	int			vSize		= velTop->mesh->nVertsTotal - velTop->bcs->size[1];
	int			eSize		= etaInt->mesh->nVertsTotal - etaInt->bcs->size[0];
	int			nksp;

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
	KSPGetPC( ksp, &pc );
	PCSetType( pc, PCFIELDSPLIT );
	PCFieldSplitSetType( pc, PC_COMPOSITE_MULTIPLICATIVE );
	ISCreateStride( MPI_COMM_WORLD, uSize, 0, 1, &utis );
	ISCreateStride( MPI_COMM_WORLD, vSize, uSize, 1, &vtis );
	ISCreateStride( MPI_COMM_WORLD, uSize, uSize + vSize, 1, &ubis );
	ISCreateStride( MPI_COMM_WORLD, vSize, uSize + vSize + uSize, 1, &vbis );
	ISCreateStride( MPI_COMM_WORLD, eSize, 2*( uSize + vSize ), 1, &eis );
	PCFieldSplitSetIS( pc, utis );
	PCFieldSplitSetIS( pc, vtis );
	PCFieldSplitSetIS( pc, ubis );
	PCFieldSplitSetIS( pc, vbis );
	PCFieldSplitSetIS( pc, eis );
	PCFieldSplitGetSubKSP( pc, &nksp, &subksp );
	MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
	if( unull ) { 
		KSPSetNullSpace( subksp[0], null ); 
		KSPSetNullSpace( subksp[2], null );
	}
	if( vnull ) { 
		KSPSetNullSpace( subksp[1], null ); 
		KSPSetNullSpace( subksp[3], null );
	}
	if( enull ) { 
		KSPSetNullSpace( subksp[4], null ); 
	}

	KSPSetOptionsPrefix( ksp, "sw2l_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );

	KSPDestroy( ksp );
	ISDestroy( utis );
	ISDestroy( vtis );
	ISDestroy( ubis );
	ISDestroy( vbis );
	ISDestroy( eis );
	MatNullSpaceDestroy( null );
}
