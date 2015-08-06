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
#include "BarotropicEqn.h"
#include "SW2L_Linear.h"

using namespace std;
using std::string;

SW2L_Linear::SW2L_Linear( Field* _velTop, Field* _velBot, Field* _etaInt, Field* _velTopPrev, Field* _velBotPrev, Field* _etaIntPrev, 
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
	MatSeqAIJSetPreallocation( A, 4*2*velTop->mesh->el->nNodes + 4*etaInt->mesh->el->nNodes, PETSC_NULL );

	VecCreate( MPI_COMM_WORLD, &b );
	VecSetSizes( b, size, PETSC_DETERMINE );
	VecSetType( b, VECSEQ );
	VecSetOption( b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );

	unull = vnull = enull = false;
}

SW2L_Linear::~SW2L_Linear() {
	delete velTopTemp;
	delete velBotTemp;
	delete etaIntTemp;
	MatDestroy( A );
	VecDestroy( b );
}

void SW2L_Linear::Assemble( int order ) {
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

	/* ei-u2 block */
	VecDuplicate( b, &x );
	row = new Vector( "row", etaInt, x, NULL, 0 );
	col = new Vector( "col", velBot, x, NULL, 0 );
	rhs = new Vector( "rhs", etaInt, b, NULL, 0 );
	diveMat = new Divergence( "ei-u2-dive-mat", etaInt, velBot, dt*bp->H );
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
	massMat = new MassMatrix( "ei-u2-mass-mat", etaInt, etaInt, a );
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

void SW2L_Linear::Solve( int order ) {
	velTopTemp->Copy( velTop );
	velBotTemp->Copy( velBot );
	etaIntTemp->Copy( etaInt );

	if( order == 1 ) { FirstOrder(); }
	if( order == 2 ) { SecondOrder(); }

	velTopPrev->Copy( velTopTemp );
	velBotPrev->Copy( velBotTemp );
	etaIntPrev->Copy( etaIntTemp );
}

void SW2L_Linear::FirstOrder() {
	Field*			negEta		= new Field( "neg-eta", etaInt->mesh, 1,  NULL );
	FieldRHS*		currRhs;
	GradFieldRHS*		gPresRhs;
	WindStressRHS*		windRhs;
	RHSOp**			ops;
	Vec			x, f;
	Vector			*utSol, *ubSol, *eiSol, *rhs;

	VecDuplicate( b, &x );
	VecDuplicate( b, &f );
	VecZeroEntries( x );
	VecZeroEntries( f );

	utSol = new Vector( "ut-sol", velTop, x, NULL, 0 );
	ubSol = new Vector( "ub-sol", velBot, x, NULL, 0 );
	eiSol = new Vector( "ei-sol", etaInt, x, NULL, 0 );

	/* u1 rhs */
	currRhs    = new FieldRHS( "vel-top-adv-rhs", velTop->mesh, 1.0, velTop );
	gPresRhs   = new GradFieldRHS( "grad-pres-rhs", velTop->mesh, -1.0*dt*tp->pFac, presSurf );
	windRhs    = new WindStressRHS( "wind-stress-rhs", velTop->mesh, dt*tp->tau, negEta, tp->k, tp->H );
	ops = new RHSOp*[3];
	ops[0] = currRhs;
	ops[1] = gPresRhs;
	ops[2] = windRhs;
	rhs = new Vector( "u1-rhs", velTop, f, ops, 3 );
	rhs->Assemble();
	delete rhs;

	/* u2 rhs */
	currRhs    = new FieldRHS( "vel-bot-adv-rhs", velBot->mesh, 1.0, velBot );
	gPresRhs   = new GradFieldRHS( "grad-pres-rhs", velBot->mesh, -1.0*dt*bp->pFac, presSurf );
	ops = new RHSOp*[2];
	ops[0] = currRhs;
	ops[1] = gPresRhs;
	rhs = new Vector( "u2-rhs", velBot, f, ops, 2 );
	rhs->Assemble();
	delete rhs;

	/* ei rhs */
	currRhs    = new FieldRHS( "eta-curr-rhs", etaInt->mesh, 1.0, etaInt );
	ops = new RHSOp*[1];
	ops[0] = currRhs;
	rhs = new Vector( "ei-rhs", etaInt, f, ops, 1 );
	rhs->Assemble();
	delete rhs;

	VecAXPY( f, 1.0, b );
	SolveLinSys( f, x );
	utSol->UpdateField();
	ubSol->UpdateField();
	eiSol->UpdateField();

	delete negEta;
	delete utSol;
	delete ubSol;
	delete eiSol;
	VecDestroy( x );
	VecDestroy( f );
}

void SW2L_Linear::SecondOrder() {
	Field*			negEta		= new Field( "neg-eta", etaInt->mesh, 1,  NULL );
	FieldRHS*		currRhs;
	FieldRHS*		prevRhs;
	GradFieldRHS*		gPresRhs;
	WindStressRHS*		windRhs;
	RHSOp**			ops;
	Vec			x, f;
	Vector			*utSol, *ubSol, *eiSol, *rhs;

	VecDuplicate( b, &x );
	VecDuplicate( b, &f );
	VecZeroEntries( x );
	VecZeroEntries( f );

	utSol = new Vector( "ut-sol", velTop, x, NULL, 0 );
	ubSol = new Vector( "ub-sol", velBot, x, NULL, 0 );
	eiSol = new Vector( "ei-sol", etaInt, x, NULL, 0 );

	/* u1 rhs */
	currRhs        = new FieldRHS( "vel-top-curr-rhs", velTop->mesh, +2.0, velTop );
	prevRhs	       = new FieldRHS( "vel-top-prev-rhs", velTop->mesh, -0.5, velTopPrev );
	gPresRhs       = new GradFieldRHS( "grad-pres-rhs", velTop->mesh, -1.0*dt*tp->pFac, presSurf );
	windRhs        = new WindStressRHS( "wind-stress-rhs", velTop->mesh, dt*tp->tau, negEta, tp->k, tp->H );
	ops = new RHSOp*[4];
	ops[0] = currRhs;
	ops[1] = prevRhs;
	ops[2] = gPresRhs;
	ops[3] = windRhs;
	rhs = new Vector( "u1-rhs", velTop, f, ops, 4 );
	rhs->Assemble();
	delete rhs;

	/* u2 rhs */
	currRhs        = new FieldRHS( "vel-bot-curr-rhs", velBot->mesh, +2.0, velBot );
	prevRhs        = new FieldRHS( "vel-bot-prev-rhs", velBot->mesh, -0.5, velBotPrev );
	gPresRhs       = new GradFieldRHS( "grad-pres-rhs", velBot->mesh, -1.0*dt*bp->pFac, presSurf );
	ops = new RHSOp*[2];
	ops[0] = currRhs;
	ops[1] = prevRhs;
	ops[2] = gPresRhs;
	rhs = new Vector( "u2-rhs", velBot, f, ops, 3 );
	rhs->Assemble();
	delete rhs;

	/* ei rhs */
	currRhs        = new FieldRHS( "eta-curr-rhs", etaInt->mesh, +2.0, etaInt );
	prevRhs        = new FieldRHS( "eta-prev-rhs", etaInt->mesh, -0.5, etaIntPrev );
	ops = new RHSOp*[2];
	ops[0] = currRhs;
	ops[1] = prevRhs;
	rhs = new Vector( "ei-rhs", etaInt, f, ops, 2 );
	rhs->Assemble();
	delete rhs;

	VecAXPY( f, 1.0, b );
	SolveLinSys( f, x );
	utSol->UpdateField();
	ubSol->UpdateField();
	eiSol->UpdateField();

	delete negEta;
	delete utSol;
	delete ubSol;
	delete eiSol;
	VecDestroy( x );
	VecDestroy( f );
}

void SW2L_Linear::SolveLinSys( Vec f, Vec x ) {
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

GLayerVector_Linear::GLayerVector_Linear( string _name, Mesh* _mesh, double _constant, Field* _field, 
			    Field* _height1, Field* _height2, Field* _velocity1, Field* _velocity2, 
			    double _H1, double _H2, double _nu, double _gamma, double _gPrime, double _tau0, double _kws, Field* _eta ) : 
RHSOp( _name, _mesh, _constant, _field ) 
{
	height1		= _height1;
	height2		= _height2;
	velocity1	= _velocity1;
	velocity2	= _velocity2;
	H1		= _H1;
	H2		= _H2;
	nu		= _nu;
	gamma		= _gamma;
	gPrime		= _gPrime;
	tau0		= _tau0;
	kws		= _kws;
	eta		= _eta;

	db     = new double*[1];  db[0]  = new double[2];
	dh1    = new double*[1];  dh1[0] = new double[2];
	dh2    = new double*[1];  dh2[0] = new double[2];
	du1    = new double*[2];  du1[0] = new double[2];  du1[1] = new double[2];
	du2    = new double*[2];  du2[0] = new double[2];  du2[1] = new double[2];
	de     = new double*[1];  de[0]  = new double[2];

	db[0][0] = db[0][1] = 0.0;
}

GLayerVector_Linear::~GLayerVector_Linear() {
	delete[] db[0];   delete[] db;
	delete[] dh1[0];  delete[] dh1;
	delete[] dh2[0];  delete[] dh2;
	delete[] du1[0];  delete[] du1[1];  delete[] du1;
	delete[] du2[0];  delete[] du2[1];  delete[] du2;
	delete[] de[0];   delete[] de;
}

void GLayerVector_Linear::AssembleElement( int el_i, double* G ) {
	double *coord, weight, detJac, *Ni, **GNix, h1, u1[2], h2, u2[2], a, b = 0.0, gCoord[2], e;

	for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
		coord  = mesh->el->quadPts[pt_i]->coord;
		weight = mesh->el->quadPts[pt_i]->weight;
		Ni     = mesh->el->ShapeFuncs( pt_i );
		GNix   = mesh->ShapeFuncDerivs( el_i, pt_i, &detJac );

		mesh->LocalToGlobal( coord, el_i, gCoord );

		if( field != NULL ) {
			field->InterpLocal( el_i, coord, &b );
			field->InterpDerivsGlobal( gCoord, db );
		}
		height1->InterpLocal( el_i, coord, &h1 );
		height2->InterpLocal( el_i, coord, &h2 );
		height1->InterpDerivsGlobal( gCoord, dh1 );
		height2->InterpDerivsGlobal( gCoord, dh2 );
		velocity1->InterpLocal( el_i, coord, u1 );
		velocity2->InterpLocal( el_i, coord, u2 );
		velocity1->InterpDerivsGlobal( gCoord, du1 );
		velocity2->InterpDerivsGlobal( gCoord, du2 );
		eta->InterpDerivsGlobal( gCoord, de );
		eta->InterpLocal( el_i, coord, &e );

		a = detJac*weight*constant/( H1 + H2 - b );

		for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
			/* viscous terms */
			G[node_i*2+0] += a*nu*H1*( 2.0*GNix[0][node_i]*du1[0][0] + GNix[1][node_i]*( du1[0][1] + du1[1][0] ) );
			G[node_i*2+1] += a*nu*H1*( 2.0*GNix[1][node_i]*du1[1][1] + GNix[0][node_i]*( du1[0][1] + du1[1][0] ) );
			G[node_i*2+0] += a*nu*H2*( 2.0*GNix[0][node_i]*du2[0][0] + GNix[1][node_i]*( du2[0][1] + du2[1][0] ) );
			G[node_i*2+1] += a*nu*H2*( 2.0*GNix[1][node_i]*du2[1][1] + GNix[0][node_i]*( du2[0][1] + du2[1][0] ) );

			/* bottom friction */
			G[node_i*2+0] += a*gamma*H2*u2[0]*Ni[node_i];
			G[node_i*2+1] += a*gamma*H2*u2[1]*Ni[node_i];

			/* bottom layer interface and topography gradient */
			G[node_i*2+0] += a*gPrime*H2*dh2[0][0]*Ni[node_i];
			G[node_i*2+1] += a*gPrime*H2*dh2[0][1]*Ni[node_i];

			/* top layer wind stress */
			G[node_i*2+0] -= a*tau0*cos( kws*gCoord[1] )*Ni[node_i];
		}
	}
}

void AssembleGFirstOrder_Linear( Params* tp, Params* bp, BarotropicEqn* bt, Field* velTop, Field* velBot, Field* etaInt, Vec* G ) {
	Vector*		rhs;
	RHSOp**		ops;
	Field*		hTop;
	Field*		hBot;
	GLayerVector_Linear*	gCurrRHS;

	CreateVector( velTop, G );
	bt->CalcHeights( tp, bp, etaInt, &hTop, &hBot );

	gCurrRHS = new GLayerVector_Linear( "g-curr", velTop->mesh, -bt->dt, NULL, hTop, hBot, velTop, velBot, 
				            tp->H, bp->H, bp->nu, bp->gamma, bp->hFac, tp->tau, tp->k, etaInt );
	ops = new RHSOp*[1];
	ops[0] = gCurrRHS;
	rhs = new Vector( "G", velTop, *G, ops, 1 );
	rhs->Assemble();

	delete rhs;
	delete hTop;
	delete hBot;
}

void AssembleGSecondOrder_Linear( Params* tp, Params* bp, BarotropicEqn* bt, Field* velTop, Field* velBot, Field* etaInt, 
				Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev, Vec* G ) 
{
	Vector*		rhs;
	RHSOp**		ops;
	Field*		hTop;
	Field*		hBot;
	Field*		hTopPrev;
	Field*		hBotPrev;
	GLayerVector_Linear*	gCurrRHS;
	GLayerVector_Linear*	gPrevRHS;

	CreateVector( velTop, G );

	bt->CalcHeights( tp, bp, etaInt, &hTop, &hBot );
	bt->CalcHeights( tp, bp, etaIntPrev, &hTopPrev, &hBotPrev );

	gCurrRHS = new GLayerVector_Linear( "g-curr", velTop->mesh, -2.0*bt->dt, NULL, hTop, hBot, velTop, velBot, 
				     tp->H, bp->H, bp->nu, bp->gamma, bp->hFac, tp->tau, tp->k, etaInt );
	gPrevRHS = new GLayerVector_Linear( "g-prev", velTop->mesh, bt->dt, NULL, hTopPrev, hBotPrev, velTopPrev, velBotPrev, 
				     tp->H, bp->H, bp->nu, bp->gamma, bp->hFac, tp->tau, tp->k, etaIntPrev );
	ops = new RHSOp*[2];
	ops[0] = gCurrRHS;
	ops[1] = gPrevRHS;
	rhs = new Vector( "G", velTop, *G, ops, 2 );
	rhs->Assemble();

	delete rhs;
	delete hTop;
	delete hBot;
	delete hTopPrev;
	delete hBotPrev;
}

void CalcVelBar_Linear( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* velBar ) {
	double *v1, *v2, hinv = 1.0/( tp->H + bp->H );

	for( int node_i = 0; node_i < velBar->mesh->nVertsTotal; node_i++ ) {
		v1 = velTop->vals[node_i];
		v2 = velBot->vals[node_i];
		velBar->vals[node_i][0] = hinv*( tp->H*v1[0] + bp->H*v2[0] );
		velBar->vals[node_i][1] = hinv*( tp->H*v1[1] + bp->H*v2[1] );
	}
}
