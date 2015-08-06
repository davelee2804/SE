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
#include "TwoLayerShallowWaterEqn.h"

using namespace std;
using std::string;

TwoLayerShallowWaterEqn::TwoLayerShallowWaterEqn( Field* _velocity1, Field* _pressure1, Field* _velocity2, Field* _pressure2, 
						  double _dt, int _svvCutoff, Field* _topography, SWParams* params, double _gPrime ) 
{
	velocity1  = _velocity1;
	pressure1  = _pressure1;
	velocity2  = _velocity2;
	pressure2  = _pressure2;

	L          = params->L;
	H          = params->H;
	U          = params->U;
	f0         = params->f0;
	beta       = params->beta;
	g          = params->g;
	gamma      = params->gamma;
	tau        = params->tau;
	kws        = params->kws;
	nu         = params->nu;
	gPrime     = _gPrime;
	dt         = _dt;
	svvCutoff  = _svvCutoff;
	topography = _topography;

	uNullSp = vNullSp = pNullSp = false;
	nNonLinIts = 1;

	InitObjs();

	presVecOffset1 = pressure1->bcs->vecOffset;
	presVecOffset2 = pressure2->bcs->vecOffset; 
}

TwoLayerShallowWaterEqn::~TwoLayerShallowWaterEqn() {
	MatDestroy( A );
	VecDestroy( b );
	VecDestroy( f );
	VecDestroy( x );
}

void TwoLayerShallowWaterEqn::Solve( Field* velPrev1, Field* presPrev1, Field* velPrev2, Field* presPrev2, int assMatOrder ) {
	if( assMatOrder == 1 ) { Assemble( 1.0 ); }
	if( assMatOrder == 2 ) { Assemble( 2.0/3.0 ); }

	_Solve( velPrev1, presPrev1, velPrev2, presPrev2 );
}

void TwoLayerShallowWaterEqn::InitObjs() {
	int vSize = 2*velocity1->mesh->nVertsTotal - velocity1->bcs->size[0] - velocity1->bcs->size[1];
	int pSize = pressure1->mesh->nVertsTotal - pressure1->bcs->size[0];
	int size  = 2*vSize + 2*pSize;
	int alloc = 4*4*velocity1->mesh->el->nNodes;

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

void TwoLayerShallowWaterEqn::Assemble( double a ) {
	/* top layer */
	MassMatrix*	v1v1MassMat;
	BetaMatrix*	v1v1BetaMat;
	Laplacian*	v1v1ViscMat;
	SVV*		v1v1SVVMat;
	Gradient*	v1p1GradMat;
	Gradient*	v1p2GradMat;
	MassMatrix*	p1p1MassMat;
	Vector*		v2SolVec;
	Vector*		p2SolVec;
	Vector*		v2RHSVec;
	Vector*		p2RHSVec;
	Matrix*		Av1v1;
	Matrix*		Av1p1;
	Matrix*		Av1p2;
	Matrix*		Ap1p1;
	/* bottom layer */
	MassMatrix*	v2v2MassMat;
	BetaMatrix*	v2v2BetaMat;
	Laplacian*	v2v2ViscMat;
	SVV*		v2v2SVVMat;
	Gradient*	v2p1GradMat;
	Gradient*	v2p2GradMat;
	MassMatrix*	p2p2MassMat;
	Vector*		v1SolVec;
	Vector*		p1SolVec;
	Vector*		v1RHSVec;
	Vector*		p1RHSVec;
	Matrix*		Av2v2;
	Matrix*		Av2p1;
	Matrix*		Av2p2;
	Matrix*		Ap2p2;
	Operator**	ops;

	v1SolVec = new Vector( "vSol1", velocity1, x, NULL, 0 );
	p1SolVec = new Vector( "pSol1", pressure1, x, NULL, 0 );
	v1RHSVec = new Vector( "vRHS1", velocity1, b, NULL, 0 );
	p1RHSVec = new Vector( "pRHS1", pressure1, b, NULL, 0 );

	v2SolVec = new Vector( "vSol2", velocity2, x, NULL, 0 );
	p2SolVec = new Vector( "pSol2", pressure2, x, NULL, 0 );
	v2RHSVec = new Vector( "vRHS2", velocity2, b, NULL, 0 );
	p2RHSVec = new Vector( "pRHS2", pressure2, b, NULL, 0 );

	/* top layer */
	v1v1MassMat = new MassMatrix( "v1Mass", velocity1, velocity1, 1.0 );
	v1v1BetaMat = new BetaMatrix( "v1Beta", velocity1, velocity1, a*dt, f0*L/U, beta*L*L/U );
	v1v1ViscMat = new Laplacian( "v1Visc", velocity1, velocity1, a*dt*nu/(L*U) );
	v1v1SVVMat  = new SVV( "v1SVV", velocity1, velocity1, a*dt/velocity1->mesh->el->N, svvCutoff );
	ops = new Operator*[4];
	ops[0] = v1v1MassMat;
	ops[1] = v1v1BetaMat;
	ops[2] = v1v1ViscMat;
	ops[3] = v1v1SVVMat;
	Av1v1 = new Matrix( "Av1v1", A, v1SolVec, v1SolVec, v1RHSVec, ops, 4 );

	v1p1GradMat = new Gradient( "v1p1Grad", velocity1, pressure1, a*dt*g*H/(U*U) );
	ops = new Operator*[1];
	ops[0] = v1p1GradMat;
	Av1p1 = new Matrix( "Av1p1", A, v1SolVec, p1SolVec, v1RHSVec, ops, 1 );

	v1p2GradMat = new Gradient( "v1p2Grad", velocity1, pressure2, a*dt*g*H/(U*U) );
	ops = new Operator*[1];
	ops[0] = v1p2GradMat;
	Av1p2 = new Matrix( "Av1p2", A, v1SolVec, p2SolVec, v1RHSVec, ops, 1 );

	p1p1MassMat = new MassMatrix( "p1p1Mass", pressure1, pressure1, 1.0 );
	ops = new Operator*[1];
	ops[0] = p1p1MassMat;
	Ap1p1 = new Matrix( "Ap1p1", A, p1SolVec, p1SolVec, p1RHSVec, ops, 1 );

	/* bottom layer */
	v2v2MassMat = new MassMatrix( "v2Mass", velocity2, velocity2, 1.0 + a*dt*gamma*L/U );
	v2v2BetaMat = new BetaMatrix( "v2Beta", velocity2, velocity2, a*dt, f0*L/U, beta*L*L/U );
	v2v2ViscMat = new Laplacian( "v2Visc", velocity2, velocity2, a*dt*nu/(L*U) );
	v2v2SVVMat  = new SVV( "v2SVV", velocity2, velocity2, a*dt/velocity1->mesh->el->N, svvCutoff );
	ops = new Operator*[4];
	ops[0] = v2v2MassMat;
	ops[1] = v2v2BetaMat;
	ops[2] = v2v2ViscMat;
	ops[3] = v2v2SVVMat;
	Av2v2 = new Matrix( "Av2v2", A, v2SolVec, v2SolVec, v2RHSVec, ops, 4 );

	v2p1GradMat = new Gradient( "v2p1Grad", velocity2, pressure1, a*dt*g*H/(U*U) );
	ops = new Operator*[1];
	ops[0] = v1p1GradMat;
	Av2p1 = new Matrix( "Av2p1", A, v2SolVec, p1SolVec, v2RHSVec, ops, 1 );

	v2p2GradMat = new Gradient( "v2p2Grad", velocity2, pressure2, a*dt*(g + gPrime)*H/(U*U) );
	ops = new Operator*[1];
	ops[0] = v2p2GradMat;
	Av2p2 = new Matrix( "Av2p2", A, v2SolVec, p2SolVec, v2RHSVec, ops, 1 );

	p2p2MassMat = new MassMatrix( "p2p2Mass", pressure2, pressure2, 1.0 );
	ops = new Operator*[1];
	ops[0] = p2p2MassMat;
	Ap2p2 = new Matrix( "Ap2p2", A, p2SolVec, p2SolVec, p2RHSVec, ops, 1 );

	VecZeroEntries( b );
	MatZeroEntries( A );

	Av1v1->Assemble();
	Av1p1->Assemble();
	Av1p2->Assemble();
	Ap1p1->Assemble();
	Av2v2->Assemble();
	Av2p1->Assemble();
	Av2p2->Assemble();
	Ap2p2->Assemble();
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );

	delete v1SolVec;
	delete p1SolVec;
	delete p1RHSVec;
	delete v1RHSVec;
	delete v2SolVec;
	delete p2SolVec;
	delete p2RHSVec;
	delete v2RHSVec;
	delete Av1v1;
	delete Av1p1;
	delete Av1p2;
	delete Ap1p1;
	delete Av2v2;
	delete Av2p1;
	delete Av2p2;
	delete Ap2p2;
}

void TwoLayerShallowWaterEqn::_Solve( Field* velPrev1, Field* presPrev1, Field* velPrev2, Field* presPrev2 ) {
	Vector*			v1SolVec;
	Vector*			p1SolVec;
	Vector*			v1RHSVec;
	Vector*			p1RHSVec;
	Vector*			v2SolVec;
	Vector*			p2SolVec;
	Vector*			v2RHSVec;
	Vector*			p2RHSVec;
	DivPhiVelMatrix*  	dpvMat1;
	OIFS*			p1NonLin;
	DivPhiVelMatrix*  	dpvMat2;
	OIFS*			p2NonLin;
	FieldRHS*		v1SLRHS;
	FieldRHS*		v2SLRHS;
	GradFieldRHS*		v1TopoRHS 	= NULL;
	GradFieldRHS*		v2TopoRHS 	= NULL;
	WindStressRHS*		v1WindRHS	= NULL;
	FieldRHS*		p1RHS;
	FieldRHS*		p2RHS;
	RHSOp**			rhs;
	KSP			ksp, *subksp;
	PC			pc;
	int			nksp;
	IS			u1is, v1is, p1is, u2is, v2is, p2is;
	MatNullSpace		null;
	int			uSize		= velocity1->mesh->nVertsTotal - velocity1->bcs->size[0];
	int			vSize		= velocity1->mesh->nVertsTotal - velocity1->bcs->size[1];
	int			pSize 		= pressure1->mesh->nVertsTotal - pressure1->bcs->size[0];
	Advector*		v1Adv;
	Advector*		v2Adv;
	double			a		= ( velPrev1 ) ? 2.0/3.0 : 1.0;
	int			nVelOps;
	int			velOp_i;
	Field*			pWindStress	= new Field( "pWindStress", pressure1->mesh, 1, NULL );

	v1SolVec = new Vector( "v1Sol", velocity1, x, NULL, 0 );
	p1SolVec = new Vector( "p1Sol", pressure1, x, NULL, 0 );
	v2SolVec = new Vector( "v2Sol", velocity2, x, NULL, 0 );
	p2SolVec = new Vector( "p2Sol", pressure2, x, NULL, 0 );

	if( velPrev1 ) { 
		v1Adv = new Advector( velocity1, velocity1, velPrev1, velPrev1 ); 
		v2Adv = new Advector( velocity2, velocity2, velPrev2, velPrev1 ); 
	}
	else { 
		v1Adv = new Advector( velocity1, velocity1 ); 
		v2Adv = new Advector( velocity2, velocity2 ); 
	}
	v1Adv->Advect( dt );
	v2Adv->Advect( dt );

	pressure1->bcs->vecOffset = 0;
	for( int node_i = 0; node_i < pressure1->mesh->nVertsTotal; node_i++ ) {
		if( pressure1->bcs->fieldToVecMap[node_i] != -1 ) {
			pressure1->bcs->fieldToVecMap[node_i] -= presVecOffset1;
		}
	}
	dpvMat1 = new DivPhiVelMatrix( "dpvMat1", pressure1, pressure1, 1.0, velocity1 );
	p1NonLin = new OIFS( dpvMat1, pressure1, presPrev1, velocity1, velPrev1, dt, nNonLinIts );
	p1NonLin->Solve();
	pressure1->bcs->vecOffset = presVecOffset1;
	for( int node_i = 0; node_i < pressure1->mesh->nVertsTotal; node_i++ ) {
		if( pressure1->bcs->fieldToVecMap[node_i] != -1 ) {
			pressure1->bcs->fieldToVecMap[node_i] += presVecOffset1;
		}
	}

	pressure2->bcs->vecOffset = 0;
	for( int node_i = 0; node_i < pressure2->mesh->nVertsTotal; node_i++ ) {
		if( pressure2->bcs->fieldToVecMap[node_i] != -1 ) {
			pressure2->bcs->fieldToVecMap[node_i] -= presVecOffset2;
		}
	}
	dpvMat2 = new DivPhiVelMatrix( "dpvMat2", pressure2, pressure2, 1.0, velocity2 );
	p2NonLin = new OIFS( dpvMat2, pressure2, presPrev2, velocity2, velPrev2, dt, nNonLinIts );
	p2NonLin->Solve();
	pressure2->bcs->vecOffset = presVecOffset2;
	for( int node_i = 0; node_i < pressure2->mesh->nVertsTotal; node_i++ ) {
		if( pressure2->bcs->fieldToVecMap[node_i] != -1 ) {
			pressure2->bcs->fieldToVecMap[node_i] += presVecOffset2;
		}
	}

	/* top layer */
	nVelOps = 1;
	velOp_i = 0;
	if( topography   ) { nVelOps++; }
	if( tau > 1.0e-6 ) { nVelOps++; }
	rhs = new RHSOp*[nVelOps];
	v1SLRHS = new FieldRHS( "v1RHS-SL", velocity1->mesh, a, v1Adv->fieldSL );
	rhs[velOp_i] = v1SLRHS;
	if( topography ) {
		velOp_i++;
		v1TopoRHS = new GradFieldRHS( "v1TopoRHS", velocity1->mesh, -a*dt*g*H/(U*U), topography );
		rhs[velOp_i] = v1TopoRHS;
	}
	if( tau > 1.0e-6 ) {
		//if( presPrev ) {
		//	for( int node_i = 0; node_i < pressure->mesh->nVertsTotal; node_i++ ) {
		//		pWindStress->vals[node_i][0] = 2.0*pressure->vals[node_i][0] - presPrev->vals[node_i][0] - 1.0;
		//	}
		//}
		//else {
			for( int node_i = 0; node_i < pressure1->mesh->nVertsTotal; node_i++ ) {
				pWindStress->vals[node_i][0] = pressure1->vals[node_i][0] - 1.0;
			}
		//}
		velOp_i++;
		v1WindRHS = new WindStressRHS( "v1WindRHS", velocity1->mesh, a*dt*tau*L/(rho*H*U*U), pWindStress, kws, 1.0 );
		rhs[velOp_i] = v1WindRHS;
	}
	v1RHSVec = new Vector( "v1RHS", velocity1, f, rhs, nVelOps );

	p1RHS = new FieldRHS( "p1RHS", pressure1->mesh, 1.0, p1NonLin->phiTilde );
	rhs = new RHSOp*[1];
	rhs[0] = p1RHS;
	p1RHSVec = new Vector( "p1RHS", pressure1, f, rhs, 1 );

	/* bottom layer */
	nVelOps = 1;
	velOp_i = 0;
	if( topography   ) { nVelOps++; }
	rhs = new RHSOp*[nVelOps];
	v2SLRHS = new FieldRHS( "v2RHS-SL", velocity2->mesh, a, v2Adv->fieldSL );
	rhs[velOp_i] = v2SLRHS;
	if( topography ) {
		velOp_i++;
		v2TopoRHS = new GradFieldRHS( "v2TopoRHS", velocity2->mesh, -a*dt*(g + gPrime)*H/(U*U), topography );
		rhs[velOp_i] = v2TopoRHS;
	}
	v2RHSVec = new Vector( "v2RHS", velocity2, f, rhs, nVelOps );

	p2RHS = new FieldRHS( "p2RHS", pressure2->mesh, 1.0, p2NonLin->phiTilde );
	rhs = new RHSOp*[1];
	rhs[0] = p2RHS;
	p2RHSVec = new Vector( "p2RHS", pressure2, f, rhs, 1 );

	VecZeroEntries( f );
	v1RHSVec->Assemble();
	p1RHSVec->Assemble();
	v2RHSVec->Assemble();
	p2RHSVec->Assemble();
	VecAXPY( f, 1.0, b );

	KSPCreate( MPI_COMM_WORLD, &ksp );
	KSPSetOperators( ksp, A, A, SAME_NONZERO_PATTERN );
	KSPGetPC( ksp, &pc );
	PCSetType( pc, PCFIELDSPLIT );
	PCFieldSplitSetType( pc, PC_COMPOSITE_MULTIPLICATIVE );
	ISCreateStride( MPI_COMM_WORLD, uSize, 0, 1, &u1is );
	ISCreateStride( MPI_COMM_WORLD, vSize, uSize, 1, &v1is );
	ISCreateStride( MPI_COMM_WORLD, pSize, uSize + vSize, 1, &p1is );
	ISCreateStride( MPI_COMM_WORLD, uSize, uSize + vSize + pSize, 1, &u2is );
	ISCreateStride( MPI_COMM_WORLD, vSize, uSize + vSize + pSize + uSize, 1, &v2is );
	ISCreateStride( MPI_COMM_WORLD, pSize, uSize + vSize + pSize + uSize + vSize, 1, &p2is );
	PCFieldSplitSetIS( pc, u1is );
	PCFieldSplitSetIS( pc, v1is );
	PCFieldSplitSetIS( pc, p1is );
	PCFieldSplitSetIS( pc, u2is );
	PCFieldSplitSetIS( pc, v2is );
	PCFieldSplitSetIS( pc, p2is );
	PCFieldSplitGetSubKSP( pc, &nksp, &subksp );
	MatNullSpaceCreate( MPI_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &null );
	if( uNullSp ) { 
		KSPSetNullSpace( subksp[0], null ); 
		KSPSetNullSpace( subksp[3], null );
	}
	if( vNullSp ) { 
		KSPSetNullSpace( subksp[1], null );
		KSPSetNullSpace( subksp[4], null );
	}
	if( pNullSp ) { 
		KSPSetNullSpace( subksp[2], null );
		KSPSetNullSpace( subksp[5], null );
	}

	KSPSetFromOptions( ksp );
	KSPSolve( ksp, f, x );
	v1SolVec->UpdateField();
	p1SolVec->UpdateField();
	v2SolVec->UpdateField();
	p2SolVec->UpdateField();

	ISDestroy( u1is );
	ISDestroy( v1is );
	ISDestroy( p1is );
	ISDestroy( u2is );
	ISDestroy( v2is );
	ISDestroy( p2is );
	MatNullSpaceDestroy( null );
	KSPDestroy( ksp );

	delete v1SolVec;
	delete p1SolVec;
	delete v1RHSVec;
	delete p1RHSVec;
	delete v2SolVec;
	delete p2SolVec;
	delete v2RHSVec;
	delete p2RHSVec;
	delete p1NonLin;
	delete p2NonLin;
	delete v1Adv;
	delete v2Adv;
	delete pWindStress;
}

void TwoLayerShallowWaterEqn::CalcEnergetics( int timeStep, ShallowWaterEnergetics* topEnergetics, ShallowWaterEnergetics* bottomEnergetics ) {
	double detJac, weight, *coord, v[2], p, gCoord[2], **d2v, h1Avg, h2Avg;
	ofstream file;
	char filename[40];

	topEnergetics->KE[timeStep] = 0.0;
	topEnergetics->PE[timeStep] = 0.0;
	topEnergetics->visc[timeStep] = 0.0;
	topEnergetics->wind[timeStep] = 0.0;
	bottomEnergetics->KE[timeStep] = 0.0;
	bottomEnergetics->PE[timeStep] = 0.0;
	bottomEnergetics->visc[timeStep] = 0.0;
	bottomEnergetics->fric[timeStep] = 0.0;

	h1Avg = pressure1->Integrate( 0, true );
	h2Avg = pressure2->Integrate( 0, true );

	d2v = new double*[2];
	d2v[0] = new double[2];
	d2v[1] = new double[2];

	/* TODO: the pressure values here assume a average layer height of 1, need to refine based on the actual height of each layer. */
	for( int el_i = 0; el_i < velocity1->mesh->nElsTotal; el_i++ ) {
		for( int pt_i = 0; pt_i < velocity1->mesh->el->nPoints; pt_i++ ) {
			coord  = velocity1->mesh->el->quadPts[pt_i]->coord;
			weight = velocity1->mesh->el->quadPts[pt_i]->weight;
			detJac = velocity1->mesh->DetJac( el_i, pt_i );
			/* top layer */
			velocity1->InterpLocal( el_i, coord, v );
			pressure1->InterpLocal( el_i, coord, &p );
			velocity1->mesh->LocalToGlobal( coord, el_i, gCoord );
			VelocitySecondDerivs( velocity1, el_i, pt_i, d2v );
			topEnergetics->KE[timeStep]   += detJac*weight*0.5*p*( v[0]*v[0] + v[1]*v[1] );
			topEnergetics->PE[timeStep]   += detJac*weight*0.5*(p - h1Avg)*(p - h1Avg);
			/* the second derivative operator is a little buggy, leaves an artifact on the element boundaries */
			topEnergetics->visc[timeStep] += detJac*weight*p*( v[0]*( d2v[0][0] + d2v[0][1] ) + v[1]*( d2v[1][0] + d2v[1][1] ) )/detJac;
			topEnergetics->wind[timeStep] += detJac*weight*v[0]*cos( M_PI*gCoord[1] );
			/* bottom layer */
			velocity2->InterpLocal( el_i, coord, v );
			pressure2->InterpLocal( el_i, coord, &p );
			velocity2->mesh->LocalToGlobal( coord, el_i, gCoord );
			VelocitySecondDerivs( velocity2, el_i, pt_i, d2v );
			bottomEnergetics->KE[timeStep]   += detJac*weight*0.5*p*( v[0]*v[0] + v[1]*v[1] );
			bottomEnergetics->PE[timeStep]   += detJac*weight*0.5*(p - h2Avg)*(p - h2Avg);
			/* the second derivative operator is a little buggy, leaves an artifact on the element boundaries */
			bottomEnergetics->visc[timeStep] += detJac*weight*p*( v[0]*( d2v[0][0] + d2v[0][1] ) + v[1]*( d2v[1][0] + d2v[1][1] ) )/detJac;
			bottomEnergetics->fric[timeStep] += detJac*weight*p*( v[0]*v[0] + v[1]*v[1] );
		}
	}
	topEnergetics->KE[timeStep]      *= rho*H*U*U;
	topEnergetics->PE[timeStep]      *= rho*g*H*H;
	topEnergetics->visc[timeStep]    *= rho*nu*H*U*U/(L*L);
	topEnergetics->wind[timeStep]    *= tau*U;
	bottomEnergetics->KE[timeStep]   *= rho*H*U*U;
	bottomEnergetics->PE[timeStep]   *= rho*g*H*H;
	bottomEnergetics->visc[timeStep] *= rho*nu*H*U*U/(L*L);
	bottomEnergetics->fric[timeStep] *= rho*gamma*H*U*U;

	delete[] d2v[0];
	delete[] d2v[1];
	delete[] d2v;

	sprintf( filename, "energetics.sw" );
	file.open( filename );
	for( int timeStep_i = 0; timeStep_i <= timeStep; timeStep_i++ ) {
		file << topEnergetics->KE[timeStep_i] << "\t" << topEnergetics->PE[timeStep_i] << "\t" << 
			bottomEnergetics->KE[timeStep_i] << "\t" << bottomEnergetics->PE[timeStep_i] << "\t" << 
			topEnergetics->visc[timeStep_i] << "\t" << bottomEnergetics->visc[timeStep_i] << 
			bottomEnergetics->fric[timeStep_i] << "\t" << topEnergetics->wind[timeStep_i] << endl;
	}
	file.close();
}
