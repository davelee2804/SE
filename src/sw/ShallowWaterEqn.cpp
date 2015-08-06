#include <iostream>
#include <fstream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SWOps.h"
#include "ShallowWaterEqn.h"

using namespace std;
using std::string;

ShallowWaterEqn::ShallowWaterEqn( Field* _velocity, Field* _pressure, double _dt, int _svvCutoff, Field* _topography, SWParams* _params ) {
	velocity   = _velocity;
	pressure   = _pressure;

	params     = _params;
	dt         = _dt;
	svvCutoff  = _svvCutoff;
	topography = _topography;

	uNullSp = vNullSp = false;
	pNullSp = true;

	qgScale    = 1.0;

	InitObjs();
}

ShallowWaterEqn::~ShallowWaterEqn() {
	MatDestroy( A );
	VecDestroy( b );
	VecDestroy( f );
	VecDestroy( x );
}

void ShallowWaterEqn::InitObjs() {
	int vSize = 2*velocity->mesh->nVertsTotal - velocity->bcs->size[0] - velocity->bcs->size[1];
	int pSize = pressure->mesh->nVertsTotal - pressure->bcs->size[0];
	int size = vSize + pSize;
	int alloc = 4*4*velocity->mesh->el->nNodes;

	MatCreate( MPI_COMM_WORLD, &A );
        MatSetSizes( A, size, size, PETSC_DETERMINE, PETSC_DETERMINE );
        MatSetType( A, MATSEQAIJ );
        MatSeqAIJSetPreallocation( A, alloc, PETSC_NULL );
        MatZeroEntries( A );

	VecCreate( MPI_COMM_WORLD, &x );
        VecSetSizes( x, size, PETSC_DETERMINE );
        VecSetType( x, VECSEQ );
        VecSetOption( x, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
        VecZeroEntries( x );

	VecDuplicate( x, &b );
	VecZeroEntries( b );

	VecDuplicate( x, &f );
	VecZeroEntries( f );
}

void ShallowWaterEqn::Assemble( int order ) {
	MassMatrix*		vMassMat;
	BetaMatrix*		vBetaMat;
	Laplacian*		vViscMat;
	SVV*			vSVVMat;
	Gradient*		vGradMat;
	Divergence*		pDiveMat;
	//DivHeightMinusTopo*	pDiveMat;
	MassMatrix*		pMassMat;
	Operator**		ops;
	Vector*			vSolVec;
	Vector*			pSolVec;
	Vector*			vRHSVec;
	Vector*			pRHSVec;
	Matrix*			Avv;
	Matrix*			Avp;
	Matrix*			Apv;
	Matrix*			App;
	double			a		= ( order == 1 ) ? 1.0 : 1.5;

	vSolVec = new Vector( "vSol", velocity, x, NULL, 0 );
	pSolVec = new Vector( "pSol", pressure, x, NULL, 0 );
	vRHSVec = new Vector( "vRHS", velocity, b, NULL, 0 );
	pRHSVec = new Vector( "pRHS", pressure, b, NULL, 0 );

	vMassMat = new MassMatrix( "vMass", velocity, velocity, a + dt*params->gamma ); /* include the frictional term here also */
	vBetaMat = new BetaMatrix( "vBeta", velocity, velocity, dt, params->f0, params->beta );
	vViscMat = new Laplacian( "vVisc", velocity, velocity, dt*params->nu );
	vSVVMat  = new SVV( "vSVV", velocity, velocity, dt/velocity->mesh->el->N, svvCutoff );
	ops = new Operator*[4];
	ops[0] = vMassMat;
	ops[1] = vBetaMat;
	ops[2] = vViscMat;
	ops[3] = vSVVMat;
	Avv = new Matrix( "Avv", A, vSolVec, vSolVec, vRHSVec, ops, 4 );

	vGradMat = new Gradient( "vGrad", velocity, pressure, dt*params->g );
	ops = new Operator*[1];
	ops[0] = vGradMat;
	Avp = new Matrix( "Avp", A, vSolVec, pSolVec, vRHSVec, ops, 1 );

	pDiveMat = new Divergence( "pDive", pressure, velocity, dt*qgScale );
	//pDiveMat = new DivHeightMinusTopo( "pDive", pressure, velocity, a*dt, params->H_i/params->H, topography );
	ops = new Operator*[1];
	ops[0] = pDiveMat;
	Apv = new Matrix( "Apv", A, pSolVec, vSolVec, pRHSVec, ops, 1 );

	pMassMat = new MassMatrix( "pMass", pressure, pressure, a );
	ops = new Operator*[1];
	ops[0] = pMassMat;
	App = new Matrix( "App", A, pSolVec, pSolVec, pRHSVec, ops, 1 );

	VecZeroEntries( b );
	MatZeroEntries( A );

	Avv->Assemble();
	Avp->Assemble();
	Apv->Assemble();
	App->Assemble();
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );

	delete vSolVec;
	delete pSolVec;
	delete pRHSVec;
	delete vRHSVec;
	delete Avv;
	delete Avp;
	delete Apv;
	delete App;
}

void ShallowWaterEqn::Solve( Field* velPrev, Field* presPrev ) {
	Field* velTemp 	= NULL;
	Field* presTemp = NULL;

	if( velPrev != NULL && presPrev != NULL ) {
		velTemp  = new Field( "vel-temp",  velocity->mesh, 2, NULL );
		presTemp = new Field( "pres-temp", pressure->mesh, 1, NULL );
		velTemp->Copy( velocity );
		presTemp->Copy( pressure );

		SecondOrder( velPrev, presPrev );

		velPrev->Copy( velTemp );
		presPrev->Copy( presTemp );
	}
	else {
		FirstOrder();
	}
}

void ShallowWaterEqn::FirstOrder() {
	Vector*			vSolVec;
	Vector*			pSolVec;
	Vector*			vRHSVec;
	Vector*			pRHSVec;
	FieldRHS*		vRHS;
	GradFieldRHS*		vTopoRHS 		= NULL;	
	WindStressRHS*		vWindRHS		= NULL;
	FieldRHS*		pRHS;
	DivHeightVelRHS*	divPresVelRHS;
	RHSOp**			rhs;
	int			nVelOps 		= 1;
	int			velOp_i			= 0;
	Advector*		vAdv			= new Advector( velocity, velocity );
	
	VecZeroEntries( x );
	vSolVec = new Vector( "vSol", velocity, x, NULL, 0 );
	pSolVec = new Vector( "pSol", pressure, x, NULL, 0 );

	/* velocity rhs */
	vAdv->Advect( dt );
	if( params->tau > 1.0e-6 ) 	{ nVelOps++; }
	if( topography != NULL )	{ nVelOps++; }	
	rhs = new RHSOp*[nVelOps];
	vRHS = new FieldRHS( "vRHS", velocity->mesh, 1.0, vAdv->fieldSL );
	rhs[velOp_i] = vRHS;
	if( topography ) {
		velOp_i++;
		vTopoRHS = new GradFieldRHS( "vTopoRHS", velocity->mesh, -dt*params->g, topography );
		rhs[velOp_i] = vTopoRHS;
	}
	if( params->tau > 1.0e-6 ) {
		velOp_i++;
		vWindRHS = new WindStressRHS( "vWindRHS", velocity->mesh, dt*params->tau, pressure, params->kws, 1.0 );
		rhs[velOp_i] = vWindRHS;
	}
	vRHSVec = new Vector( "vRHS", velocity, f, rhs, nVelOps );

	/* pressure rhs */
	pRHS = new FieldRHS( "pRHS", pressure->mesh, 1.0, pressure );
	divPresVelRHS = new DivHeightVelRHS( "div-pres-vel-rhs", pressure->mesh, -dt, pressure, velocity );
	rhs = new RHSOp*[2];
	rhs[0] = pRHS;
	rhs[1] = divPresVelRHS;
	pRHSVec = new Vector( "pRHS", pressure, f, rhs, 2 );

	VecZeroEntries( f );
	vRHSVec->Assemble();
	pRHSVec->Assemble();
	VecAXPY( f, 1.0, b );

	SolveLinAlg( f, x );

	vSolVec->UpdateField();
	pSolVec->UpdateField();

	delete vSolVec;
	delete pSolVec;
	delete vRHSVec;
	delete pRHSVec;
	delete vAdv;
}

void ShallowWaterEqn::SecondOrder( Field* velPrev, Field* presPrev ) {
	Vector*			vSolVec;
	Vector*			pSolVec;
	Vector*			vRHSVec;
	Vector*			pRHSVec;
	FieldRHS*		vRHS;
	GradFieldRHS*		vTopoRHS 		= NULL;	
	WindStressRHS*		vWindRHS		= NULL;
	FieldRHS*		pRHS;
	FieldRHS*		pPrevRHS;
	DivHeightVelRHS*	divPresVelRHS;
	DivHeightVelRHS*	divPresVelPrevRHS;
	RHSOp**			rhs;
	int			nVelOps 		= 1;
	int			velOp_i			= 0;
	Advector*		vAdv			= new Advector( velocity, velocity, velPrev, velPrev );
	
	VecZeroEntries( x );
	vSolVec = new Vector( "vSol", velocity, x, NULL, 0 );
	pSolVec = new Vector( "pSol", pressure, x, NULL, 0 );

	/* velocity rhs */
	vAdv->Advect( dt );
	if( params->tau > 1.0e-6 ) 	{ nVelOps++; }
	if( topography != NULL )	{ nVelOps++; }	
	rhs = new RHSOp*[nVelOps];
	vRHS = new FieldRHS( "vRHS", velocity->mesh, 1.0, vAdv->fieldSL );
	rhs[velOp_i] = vRHS;
	if( topography ) {
		velOp_i++;
		vTopoRHS = new GradFieldRHS( "vTopoRHS", velocity->mesh, -dt*params->g, topography );
		rhs[velOp_i] = vTopoRHS;
	}
	if( params->tau > 1.0e-6 ) {
		velOp_i++;
		vWindRHS = new WindStressRHS( "vWindRHS", velocity->mesh, dt*params->tau, pressure, params->kws, 1.0 );
		rhs[velOp_i] = vWindRHS;
	}
	vRHSVec = new Vector( "vRHS", velocity, f, rhs, nVelOps );

	/* pressure rhs */
	pRHS = new FieldRHS( "pRHS", pressure->mesh, 2.0, pressure );
	pPrevRHS = new FieldRHS( "pPrevRHS", pressure->mesh, -0.5, presPrev );
	divPresVelRHS = new DivHeightVelRHS( "div-pres-vel-rhs", pressure->mesh, -2.0*dt, pressure, velocity );
	divPresVelPrevRHS = new DivHeightVelRHS( "div-pres-vel-prev-rhs", pressure->mesh, dt, presPrev, velPrev );
	rhs = new RHSOp*[4];
	rhs[0] = pRHS;
	rhs[1] = pPrevRHS;
	rhs[2] = divPresVelRHS;
	rhs[3] = divPresVelPrevRHS;
	pRHSVec = new Vector( "pRHS", pressure, f, rhs, 4 );

	VecZeroEntries( f );
	vRHSVec->Assemble();
	pRHSVec->Assemble();
	VecAXPY( f, 1.0, b );

	SolveLinAlg( f, x );

	vSolVec->UpdateField();
	pSolVec->UpdateField();

	delete vSolVec;
	delete pSolVec;
	delete vRHSVec;
	delete pRHSVec;
	delete vAdv;
}

void ShallowWaterEqn::SolveLinAlg( Vec f, Vec x ) {
	KSP                     ksp, *subksp;
	PC                      pc;
	int                     nksp;
	IS                      uis, vis, pis;
	MatNullSpace            null;
	int                     uSize           = velocity->mesh->nVertsTotal - velocity->bcs->size[0];
	int                     vSize           = velocity->mesh->nVertsTotal - velocity->bcs->size[1];
	int                     pSize           = pressure->mesh->nVertsTotal - pressure->bcs->size[0];

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
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
	if( uNullSp ) { KSPSetNullSpace( subksp[0], null ); }
	if( vNullSp ) { KSPSetNullSpace( subksp[1], null ); }
	if( pNullSp ) { KSPSetNullSpace( subksp[2], null ); }

	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );

	ISDestroy( uis );
	ISDestroy( vis );
	ISDestroy( pis );
	MatNullSpaceDestroy( null );
	KSPDestroy( ksp );
}
