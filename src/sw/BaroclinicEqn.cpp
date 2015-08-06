#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "LinSys.h"
#include "SWOps.h"
#include "ShallowWaterEqn.h"
#include "BaroclinicEqn.h"

using namespace std;
using std::string;

/*
Solves the two-layer shallow water system with bottom topography subject to a rigid lid top
boundary as enforced by a surface pressure constraint

References:
	Dukowicz, J. D., R. D. Smith and R. C. Malone (1993) "A Reformulation and Implementation
	of the Bryan-Cox-Semtner Ocean Model on the Connection Machine" J. Atmos. Oceanic. Technol.
	10, 195-208
	Hallberg, R. (1997) "Stable Split Time Stepping Schemes for Large Scale Ocean Modeling"
	J. Comp. Phys. 135, 54-65
	McWilliams, J. C. (2006) Fundamentals of Geophysical Fluid Dynamics, Cambridge University
	Press
*/

BaroclinicEqn::BaroclinicEqn( Field* 			_v1Curr, 
			      Field* 			_v1Prev, 
			      Field* 			_v2Curr, 
			      Field* 			_v2Prev,
			      Field* 			_etaCurr, 
			      Field* 			_etaPrev, 
			      Field* 			_psCurr, 
			      Field* 			_psPrev,
			      ShallowWaterParams* 	_topParams, 
			      ShallowWaterParams* 	_botParams,
			      Field* 			_topo, 
			      int    			_svvCutoff,
			      double 			_dt ) 
{
	v1Curr 		= _v1Curr;
	v1Prev 		= _v1Prev;
	v2Curr 		= _v2Curr;
	v2Prev 		= _v2Prev;
	etaCurr 	= _etaCurr;
	etaPrev 	= _etaPrev;
	psCurr 		= _psCurr;
	psPrev 		= _psPrev;
	topParams 	= _topParams;
	botParams 	= _botParams;
	topo 		= _topo;
	svvCutoff 	= _svvCutoff;
	dt 		= _dt;

	/* initialise the auxillary fields */
	btVel  = new Field( "barotropic-velocity", v1Curr->mesh, 2, v1Curr->bcs );
	bcVel1 = new Field( "baroclinic-velocity-top-layer", v1Curr->mesh, 2, v1Curr->bcs );
	bcVel2 = new Field( "baroclinic-velocity-bottom-layer", v1Curr->mesh, 2, v1Curr->bcs );

	husqN = new Field( "kinetic-energy-flux-previous-time", v1Curr->mesh, 1, NULL );
	husqP = new Field( "kinetic-energy-flux-predictor-step", v1Curr->mesh, 1, NULL );

	topHeight = new Field( "top-layer-height", etaCurr->mesh, 1, NULL );
	botHeight = new Field( "bottom-layer-height", etaCurr->mesh, 1, NULL );

	BaroclinicPredictorSetup();
	AssembleBarotropic();

	swSolver = new ShallowWaterEqn( v1Curr, etaCurr, dt, v1Curr->mesh->el->N/2, topo, botParams );
	swAssembled = false;
}

BaroclinicEqn::~BaroclinicEqn() {
	delete btVel;
	delete bcVel1;
	delete bcVel2;
	delete husqN;
	delete husqP;
	delete topHeight;
	delete botHeight;

	/* free the implicit barotropic solve objects */
	MatDestroy( bt_Kinv );
	MatDestroy( bt_G );
	MatDestroy( bt_D );
	VecDestroy( bt_b );

	/* free the baroclinic predictor step objects */
	MatDestroy( bc_M_vp );
	MatDestroy( bc_M_ep );
	VecDestroy( bc_b_v );

	delete scSolver;
	delete swSolver;
	delete v1Adv;
	delete v2Adv;
}

void BaroclinicEqn::Solve( int order ) {
	Field*	bcv1p	= new Field( "bcv1p", v1Curr->mesh, 2, v1Curr->bcs );
	Field*	bcv2p	= new Field( "bcv2p", v2Curr->mesh, 2, v2Curr->bcs );
	Field*	bcep	= new Field( "bcep", etaCurr->mesh, 1, etaCurr->bcs );
	Field*	btvp 	= new Field( "btvp", btVel->mesh, 2, btVel->bcs );
	Field*	btpsp	= new Field( "btpsp", psCurr->mesh, 1, psCurr->bcs );
	Field*	v1p	= new Field( "v1p", v1Curr->mesh, 2, v1Curr->bcs );
	Field*	v2p	= new Field( "v2p", v2Curr->mesh, 2, v2Curr->bcs );
	Field*	btvNew	= new Field( "btvNew", btVel->mesh, 2, btVel->bcs );
	Field*	psNew 	= new Field( "psNew", psCurr->mesh, 1, psCurr->bcs );

	if( order == 1 ) {
		v1Adv = new Advector( v1Curr, v1Curr );
		v2Adv = new Advector( v2Curr, v2Curr );
	}
	if( order == 2 ) {
		v1Adv = new Advector( v1Curr, v1Curr, v1Prev, v1Prev );
		v2Adv = new Advector( v2Curr, v2Curr, v2Prev, v2Prev );
	}

	/* predictor step */
	CalcKineticEnergyFlux( etaCurr, v1Curr, v2Curr, husqN );
	SolveBarotropic( btvp, btpsp, husqN );
	BaroclinicPredictor( bcv1p, bcv2p, bcep );

	for( int node_i = 0; node_i < btVel->mesh->nVertsTotal; node_i++ ) {
		v1p->vals[node_i][0] = btvp->vals[node_i][0] + bcv1p->vals[node_i][0];
		v1p->vals[node_i][1] = btvp->vals[node_i][1] + bcv1p->vals[node_i][1];
		v2p->vals[node_i][0] = btvp->vals[node_i][0] + bcv2p->vals[node_i][0];
		v2p->vals[node_i][1] = btvp->vals[node_i][1] + bcv2p->vals[node_i][1];
		
	}
	CalcKineticEnergyFlux( bcep, v1p, v2p, husqP );

	/* corrector step */
	for( int node_i = 0; node_i < husqP->mesh->nVertsTotal; node_i++ ) {
		husqP->vals[node_i][0] *= 0.5;
		husqP->vals[node_i][0] += 0.5*husqN->vals[node_i][0];
	}
	for( int node_i = 0; node_i < btVel->mesh->nVertsTotal; node_i++ ) {
		btvp->vals[node_i][0] *= 0.5;
		btvp->vals[node_i][0] += 0.5*btVel->vals[node_i][0];
	}

	SolveBarotropic( btvNew, psNew, husqP );
	if( order == 1 ) {
		swSolver->Assemble( 1.0 );
	}
	if( order == 2 && !swAssembled ) {
		swSolver->Assemble( 2.0/3.0 );
		swAssembled = true;
	}
	SolveShallowWater( btvp, husqP, order );

	delete bcv1p;
	delete bcv2p;
	delete bcep;
	delete btvp;
	delete btpsp;
	delete v1p;
	delete v2p;
	delete btvNew;
	delete psNew;
}

void BaroclinicEqn::AssembleBarotropic() {
	BetaInvMatrix*		betaInvMat;
	Gradient*		gradMat;
	BarotropicDivergence*	divMat;
	Operator**		ops;
	Vec			v;
	Vec			p;
	Vector*			vSol;
	Vector*			vRhs;
	Vector*			pSol;
	Matrix*			kMat;
	Matrix*			gMat;
	Matrix*			dMat;

	CreateMatrix( btVel, btVel, &bt_Kinv );
	CreateMatrix( btVel, psCurr, &bt_G );
	CreateMatrix( psCurr, btVel, &bt_D );

	CreateVector( btVel, &bt_b );
	CreateVector( btVel, &v );
	CreateVector( psCurr, &p );

	vSol = new Vector( "velocity-solution-vector", btVel, v, NULL, 0 );
	pSol = new Vector( "pressure-solution-vector", psCurr, p, NULL, 0 );
	vRhs = new Vector( "rhs-vector", btVel, bt_b, NULL, 0 );

	/* inverse of the mass+coriolis matrix, 
	 *	    [  1    -dt.f ]			       1     [  1    +dt.f ]
	 * 	A = |             |		A^{-1} = ------------|             |
	 *	    [ +dt.f    1  ]			 1 + (dt.f)^2[ -dt.f    1  ]
	 */
	betaInvMat = new BetaInvMatrix( "beta-operator", btVel, btVel, dt, 0.0, topParams->f0, topParams->beta );
	ops = new Operator*[1];
	ops[0] = betaInvMat;
	kMat = new Matrix( "advection-coriolis-matrix", bt_Kinv, vSol, vSol, vRhs, ops, 1 );
	kMat->Assemble();

	gradMat = new Gradient( "gradient-operator", btVel, psCurr, dt/topParams->rho );
	ops = new Operator*[1];
	ops[0] = gradMat;
	gMat = new Matrix( "gradient-matrix", bt_G, vSol, pSol, vRhs, ops, 1 );	
	gMat->Assemble();

	divMat = new BarotropicDivergence( "divergence-operator", psCurr, btVel, 1.0, topo );
	ops = new Operator*[1];
	ops[0] = divMat;
	dMat = new Matrix( "divergence-matrix", bt_D, pSol, vSol, NULL, ops, 1 );
	dMat->Assemble();

	VecDestroy( v );
	VecDestroy( p );
	delete vSol;
	delete pSol;
	delete vRhs;
	delete kMat;
	delete gMat;
	delete dMat;

	scSolver = new SchurComplement( bt_Kinv, bt_G, bt_D, NULL );
}

void BaroclinicEqn::SolveBarotropic( Field* btVelP, Field* psP, Field* husq ) {
	GradFieldRHS*	dhusqRHS;
	FieldRHS*	velRHS;
	RHSOp**		rhsOps;
	Vec		v;
	Vec		p;
	Vec		f;
	Vector*		vSol;
	Vector*		pSol;
	Vector*		vRhs;

	CreateVector( btVel, &v );
	CreateVector( btVel, &f );
	CreateVector( psCurr, &p );

	cout << "solving barotropic system...\n";

	vSol = new Vector( "v", btVelP, v, NULL, 0 );
	pSol = new Vector( "p", psP, p, NULL, 0 );
	velRHS = new FieldRHS( "bt-vel-prev-rhs", btVel->mesh, 1.0, btVel );
	dhusqRHS = new GradFieldRHS( "bt-vel-dhusq-rhs", btVel->mesh, -dt, husq );
	rhsOps = new RHSOp*[2];
	rhsOps[0] = velRHS;
	rhsOps[1] = dhusqRHS;
	vRhs = new Vector( "f", btVel, f, rhsOps, 2 );
	vRhs->Assemble();

	scSolver->Solve( v, p, f, NULL );
	vSol->UpdateField();
	pSol->UpdateField();

	delete vSol;
	delete pSol;
	delete vRhs;

	VecDestroy( v );
	VecDestroy( f );
	VecDestroy( p );
}

void BaroclinicEqn::BaroclinicPredictorSetup() {
	MassMatrix*	vMassMat;
	MassMatrix*	eMassMat;
	Operator**	ops;
	Vec		v;
	Vec		e;
	Vector*		vSol;
	Vector*		eSol;
	Vector*		vRhs;
	Matrix*		vMat;
	Matrix*		eMat;

	CreateMatrix( v1Curr, v1Curr, &bc_M_vp );
	CreateMatrix( etaCurr, etaCurr, &bc_M_ep );
	CreateVector( v1Curr, &bc_b_v );

	CreateVector( v1Curr, &v );
	CreateVector( etaCurr, &e );

	vSol = new Vector( "velocity-solution-vector", v1Curr, v, NULL, 0 );
	eSol = new Vector( "interface-solution-vector", etaCurr, e, NULL, 0 );
	vRhs = new Vector( "velocity-rhs-vector", v1Curr, bc_b_v, NULL, 0 );

	vMassMat = new MassMatrix( "velocity-predictor-operator", v1Curr, v1Curr, 1.0 );
	ops = new Operator*[1];
	ops[0] = vMassMat;
	vMat = new Matrix( "velocity-predictor-matrix", bc_M_vp, vSol, vSol, vRhs, ops, 1 );
	vMat->Assemble();

	eMassMat = new MassMatrix( "interface-predictor-operator", etaCurr, etaCurr, 1.0 );
	ops = new Operator*[1];
	ops[0] = eMassMat;
	eMat = new Matrix( "interface-predictor-matrix", bc_M_ep, eSol, eSol, NULL, ops, 1 );
	eMat->Assemble();

	VecDestroy( v );
	VecDestroy( e );
	delete vSol;
	delete eSol;
	delete vRhs;
	delete vMat;
	delete eMat;
}

void BaroclinicEqn::CalcKineticEnergyFlux( Field* eta, Field* velocity1, Field* velocity2, Field* husq ) {
	double e, b, h1, h2, v1[2], v2[2], H;

	for( int node_i = 0; node_i < husq->mesh->nVertsTotal; node_i++ ) {
		eta->InterpGlobal( husq->mesh->verts[node_i], &e );
		velocity1->InterpGlobal( husq->mesh->verts[node_i], v1 );
		velocity2->InterpGlobal( husq->mesh->verts[node_i], v2 );
		topo->InterpGlobal( husq->mesh->verts[node_i], &b );
		H = topParams->H + botParams->H - b;
		h1 = topParams->H - e;
		h2 = botParams->H + e - b;
		husq->vals[node_i][0] = ( h1*( v1[0]*v1[0] + v1[1]*v1[1] ) + h2*( v2[0]*v2[0] + v2[1]*v2[1] ) )/H;
	}
}

void BaroclinicEqn::CalcHeightFields( Field* eta ) {
	double b, e;

	if( eta->mesh != topHeight->mesh || eta->mesh != botHeight->mesh ) { 
		cout << "ERROR: height field mesh not set correctly\n"; 
		exit(0);
	}

	for( int node_i = 0; node_i < eta->mesh->nVertsTotal; node_i++ ) {
		topo->InterpGlobal( eta->mesh->verts[node_i], &b );
		e = eta->vals[node_i][0];

		topHeight->vals[node_i][0] = topParams->H - e;
		botHeight->vals[node_i][0] = botParams->H + e - b;
	}
}

void BaroclinicEqn::BaroclinicPredictor( Field* bcVel1P, Field* bcVel2P, Field* etaP ) {
	double		e, t, f, tau, *coord;
	double		**dhusq, **de, **db;
	Field*		rhsFieldTop 	= new Field( "rhs-field-top", v1Curr->mesh, 2, NULL );
	Field*		rhsFieldBot	= new Field( "rhs-field-bottom", v2Curr->mesh, 2, NULL );
	Field*		rhsFieldEta	= new Field( "rhs-field-eta", etaCurr->mesh, 2, NULL );
	Vec		b;
	Vec		x;
	Vector*		rhs;
	Vector*		sol;
	FieldRHS*	fieldRHS;
	DivVelRHS*	divRHS;
	RHSOp**		rhsOp;
	KSP		ksp;
	
	dhusq    = new double*[1];
	dhusq[1] = new double[2];
	de       = new double*[1];
	de[1]    = new double[2]; 
	db       = new double*[1];
	db[1]    = new double[2]; 

	/* create the rhs fields from the current time step fields */
	for( int node_i = 0; node_i < v1Curr->mesh->nVertsTotal; node_i++ ) {
		coord = v1Curr->mesh->verts[node_i];
		etaCurr->InterpGlobal( coord, &e );
		topo->InterpGlobal( coord, &t );
		f = topParams->f0 + topParams->beta*coord[1];
		tau = topParams->tau*cos( topParams->kws*coord[1] );
		husqN->InterpDerivsGlobal( coord, dhusq );
		etaCurr->InterpDerivsGlobal( coord, de );
		topo->InterpDerivsGlobal( coord, db );

		rhsFieldTop->vals[node_i][0] = ( v1Curr->vals[node_i][0] - btVel->vals[node_i][0] ) +
					       dt*f*( v1Curr->vals[node_i][1] - btVel->vals[node_i][1] ) + 
					       v1Adv->fieldSL->vals[node_i][0] + /*??*/
					       dt*dhusq[0][0] + 
					       dt*tau/( topParams->rho*( topParams->H - e ) );

		rhsFieldTop->vals[node_i][1] = ( v1Curr->vals[node_i][1] - btVel->vals[node_i][1] ) -
					       dt*f*( v1Curr->vals[node_i][0] - btVel->vals[node_i][0] ) + 
					       v1Adv->fieldSL->vals[node_i][1] + /*??*/
					       dt*dhusq[0][1] + 
					       dt*tau/( topParams->rho*( topParams->H - e ) );

		rhsFieldBot->vals[node_i][0] = ( v2Curr->vals[node_i][0] - btVel->vals[node_i][0] ) +
					       dt*f*( v2Curr->vals[node_i][1] - btVel->vals[node_i][1] ) + 
					       v2Adv->fieldSL->vals[node_i][0] + /*??*/
					       dt*dhusq[0][0] - 
					       dt*botParams->gamma*v2Curr->vals[node_i][0] - 
					       dt*botParams->g*( de[0][0] + db[0][0] );

		rhsFieldBot->vals[node_i][1] = ( v2Curr->vals[node_i][1] - btVel->vals[node_i][1] ) -
					       dt*f*( v2Curr->vals[node_i][0] - btVel->vals[node_i][0] ) + 
					       v2Adv->fieldSL->vals[node_i][1] + /*??*/
					       dt*dhusq[0][1] - 
					       dt*botParams->gamma*v2Curr->vals[node_i][1] - 
					       dt*botParams->g*( de[0][1] + db[0][1] );

		rhsFieldEta->vals[node_i][0] = -dt*( botParams->H + e - t )*v2Curr->vals[node_i][0];
		rhsFieldEta->vals[node_i][1] = -dt*( botParams->H + e - t )*v2Curr->vals[node_i][1];
	}

	delete[] dhusq[1];
	delete[] dhusq;
	delete[] de[1];
	delete[] de;
	delete[] db[1];
	delete[] db;

	/* solve for top layer predictor step */
	CreateVector( v1Curr, &x );
	CreateVector( v1Curr, &b );

	cout << "solving for top layer velocity predictor step...\n";

	sol = new Vector( "x", bcVel1P, x, NULL, 0 );
	/* assemble the rhs vector */
	fieldRHS = new FieldRHS( "top-vel-pred-rhs", v1Curr->mesh, 1.0, rhsFieldTop );
	rhsOp = new RHSOp*[1];
	rhsOp[0] = fieldRHS;
	rhs = new Vector( "b", bcVel1P, b, rhsOp, 1 );
	rhs->Assemble();
	/* add in the boundary conditions */
	VecAXPY( b, 1.0, bc_b_v );
	/* solve */
	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, bc_M_vp, bc_M_vp, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "vel_pred_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, b, x );
	sol->UpdateField();

	/* replace the top rhs vector with the bottom rhs vector and solve for the bottom layer predictor step */
	VecZeroEntries( x );
	VecZeroEntries( b );
	delete sol;
	delete rhs;

	cout << "solving for bottom layer velocity predictor step...\n";

	sol = new Vector( "x", bcVel2P, x, NULL, 0 );
	fieldRHS = new FieldRHS( "bot-vel-pred-rhs", v2Curr->mesh, 1.0, rhsFieldBot );
	rhsOp = new RHSOp*[1];
	rhsOp[0] = fieldRHS;
	rhs = new Vector( "b", bcVel2P, b, rhsOp, 1 );
	rhs->Assemble();
	VecAXPY( b, 1.0, bc_b_v );
	KSPSolve( ksp, b, x );
	sol->UpdateField();
	
	KSPDestroy( ksp );
	VecDestroy( x );
	VecDestroy( b );
	delete sol;
	delete rhs;

	/* solve for the interface height predictor step */
	CreateVector( etaCurr, &x );
	CreateVector( etaCurr, &b );

	cout << "solving for interface predictor step...\n";

	sol = new Vector( "x", etaP, x, NULL, 0 );
	divRHS = new DivVelRHS( "eta-pred-rhs", etaCurr->mesh, 1.0, rhsFieldEta );
	rhsOp = new RHSOp*[1];
	rhsOp[0] = divRHS;
	rhs = new Vector( "b", etaP, b, rhsOp, 1 );
	rhs->Assemble();
	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, bc_M_ep, bc_M_ep, SAME_NONZERO_PATTERN );
	KSPSetOptionsPrefix( ksp, "eta_pred_" );
	KSPSetFromOptions( ksp );
	KSPSolve( ksp, b, x );
	sol->UpdateField();

	KSPDestroy( ksp );
	VecDestroy( x );
	VecDestroy( b );
	delete sol;
	delete rhs;

	delete rhsFieldTop;
	delete rhsFieldBot;
	delete rhsFieldEta;
}

/* still working on this function */
void BaroclinicEqn::SolveShallowWater( Field* btv, Field* husq, int order ) {
	Vector*			vSolVec;
	Vector*			pSolVec;
	Vector*			vRHSVec;
	Vector*			pRHSVec;
	DivPhiVelMatrix*  	dpvMat;
	OIFS*			pNonLin;
	FieldRHS*		vRHS;
	GradFieldRHS*		vTopoRHS;
	FieldRHS*		fricRHS;
	GradFieldRHS*		husqRHS;
	FieldRHS*		pRHS;
	RHSOp**			rhs;
	KSP			ksp, *subksp;
	PC			pc;
	int			nksp;
	IS			uis, vis, pis;
	MatNullSpace		null;
	int			uSize		= v2Curr->mesh->nVertsTotal - v2Curr->bcs->size[0];
	int			vSize		= v2Curr->mesh->nVertsTotal - v2Curr->bcs->size[1];
	int			pSize 		= etaCurr->mesh->nVertsTotal - etaCurr->bcs->size[0];
	double			a		= ( order == 2 ) ? 2.0/3.0 : 1.0;

	vSolVec = new Vector( "vSol", v2Curr, swSolver->x, NULL, 0 );
	pSolVec = new Vector( "pSol", etaCurr, swSolver->x, NULL, 0 );

	etaCurr->bcs->vecOffset = 0;
	for( int node_i = 0; node_i < etaCurr->mesh->nVertsTotal; node_i++ ) {
		if( etaCurr->bcs->fieldToVecMap[node_i] != -1 ) {
			etaCurr->bcs->fieldToVecMap[node_i] -= swSolver->presVecOffset;
		}
	}

	dpvMat = new DivPhiVelMatrix( "dpvMat", etaCurr, etaCurr, v2Curr );
	pNonLin = new OIFS( dpvMat, etaCurr, etaPrev, v2Curr, v2Prev, dt, 1 );
	pNonLin->Solve();

	etaCurr->bcs->vecOffset = swSolver->presVecOffset;
	for( int node_i = 0; node_i < etaCurr->mesh->nVertsTotal; node_i++ ) {
		if( etaCurr->bcs->fieldToVecMap[node_i] != -1 ) {
			etaCurr->bcs->fieldToVecMap[node_i] += swSolver->presVecOffset;
		}
	}

	/* momentum eqn rhs */
	rhs = new RHSOp*[4];
	vRHS = new FieldRHS( "vRHS", v2Curr->mesh, a, v2Adv->fieldSL );
	vTopoRHS = new GradFieldRHS( "vTopoRHS", v2Curr->mesh, -a*dt*botParams->g*botParams->H/(botParams->U*botParams->U), topo );
	fricRHS = new FieldRHS( "fricRHS", v2Curr->mesh, -a*dt*botParams->gamma, btv );
	husqRHS = new GradFieldRHS( "husqRHS", v2Curr->mesh, a*dt, husq );
	rhs[0] = vRHS;
	rhs[1] = vTopoRHS;
	rhs[2] = fricRHS;
	rhs[3] = husqRHS;
	vRHSVec = new Vector( "vRHS", v2Curr, swSolver->f, rhs, 4 );

	/* mass eqn rhs TODO */
	pRHS = new FieldRHS( "pRHS", etaCurr->mesh, 1.0, pNonLin->phiTilde );
	rhs = new RHSOp*[1];
	rhs[0] = pRHS;
	pRHSVec = new Vector( "pRHS", etaCurr, swSolver->f, rhs, 1 );

	VecZeroEntries( swSolver->f );
	vRHSVec->Assemble();
	pRHSVec->Assemble();
	VecAXPY( swSolver->f, 1.0, swSolver->b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOptionsPrefix( ksp, "shallow_water_" );
	KSPSetOperators( ksp, swSolver->A, swSolver->A, SAME_NONZERO_PATTERN );
	KSPGetPC( ksp, &pc );
	PCSetType( pc, PCFIELDSPLIT );
	PCFieldSplitSetType( pc, PC_COMPOSITE_MULTIPLICATIVE );
	ISCreateStride( MPI_COMM_WORLD, uSize, 0, 1, &uis );
	ISCreateStride( MPI_COMM_WORLD, vSize, uSize, 1, &vis );
	ISCreateStride( MPI_COMM_WORLD, pSize, uSize + vSize, 1, &pis );
	PCFieldSplitSetIS( pc, uis );
	PCFieldSplitSetIS( pc, vis );
	PCFieldSplitSetIS( pc, pis );
	PCFieldSplitGetSubKSP( pc, &nksp, &subksp );
	MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
	if( swSolver->uNullSp ) { KSPSetNullSpace( subksp[0], null ); }
	if( swSolver->vNullSp ) { KSPSetNullSpace( subksp[1], null ); }
	if( swSolver->pNullSp ) { KSPSetNullSpace( subksp[2], null ); }

	KSPSetFromOptions( ksp );
	KSPSolve( ksp, swSolver->f, swSolver->x );
	vSolVec->UpdateField();
	pSolVec->UpdateField();

	ISDestroy( uis );
	ISDestroy( vis );
	ISDestroy( pis );
	MatNullSpaceDestroy( null );
	KSPDestroy( ksp );

	delete vSolVec;
	delete pSolVec;
	delete vRHSVec;
	delete pRHSVec;
	delete pNonLin;
}
