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
#include "LinearShallowWaterEqn.h"

using namespace std;
using std::string;

LinearShallowWaterEqn::LinearShallowWaterEqn( Field* _velocity, Field* _pressure, double _f0, double _beta, double _g, 
				  double _gamma, double _tau, double _nu, double _dt, int _svvCutoff, Field* _topography ) 
{
	velocity   = _velocity;
	pressure   = _pressure;
	f0         = _f0;
	beta       = _beta;
	g          = _g;
	gamma      = _gamma;
	tau        = _tau;
	nu         = _nu;
	dt         = _dt;
	svvCutoff  = _svvCutoff;
	topography = _topography;

	uNullSp = vNullSp = pNullSp = false;

	InitObjs();
}

LinearShallowWaterEqn::~LinearShallowWaterEqn() {
	MatDestroy( A );
	VecDestroy( b );
	VecDestroy( f );
	VecDestroy( x );
}

void LinearShallowWaterEqn::Solve( Field* velPrev, Field* presPrev, int assMatOrder ) {
	if( assMatOrder == 1 ) { Assemble( 1.0 ); }
	if( assMatOrder == 2 ) { Assemble( 2.0/3.0 ); }

	_Solve( velPrev, presPrev );
}

void LinearShallowWaterEqn::InitObjs() {
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

void LinearShallowWaterEqn::Assemble( double a ) {
	MassMatrix*	vMassMat;
	BetaMatrix*	vBetaMat;
	Laplacian*	vViscMat;
	SVV*		vSVVMat;
	Gradient*	vGradMat;
	Divergence*	pDiveMat;
	MassMatrix*	pMassMat;
	Operator**	ops;
	Vector*		vSolVec;
	Vector*		pSolVec;
	Vector*		vRHSVec;
	Vector*		pRHSVec;
	Matrix*		Avv;
	Matrix*		Avp;
	Matrix*		Apv;
	Matrix*		App;

	vSolVec = new Vector( "vSol", velocity, x, NULL, 0 );
	pSolVec = new Vector( "pSol", pressure, x, NULL, 0 );
	vRHSVec = new Vector( "vRHS", velocity, b, NULL, 0 );
	pRHSVec = new Vector( "pRHS", pressure, b, NULL, 0 );

	vMassMat = new MassMatrix( "vMass", velocity, velocity, 1.0 + a*dt*gamma ); /* include the frictional term here also */
	vBetaMat = new BetaMatrix( "vBeta", velocity, velocity, a*dt, f0, beta );
	vViscMat = new Laplacian( "vVisc", velocity, velocity, a*dt*nu );
	vSVVMat  = new SVV( "vSVV", velocity, velocity, a*dt, svvCutoff );
	ops = new Operator*[4];
	ops[0] = vMassMat;
	ops[1] = vBetaMat;
	ops[2] = vViscMat;
	ops[3] = vSVVMat;
	Avv = new Matrix( "Avv", A, vSolVec, vSolVec, vRHSVec, ops, 4 );

	vGradMat = new Gradient( "vGrad", velocity, pressure, a*dt*g );
	ops = new Operator*[1];
	ops[0] = vGradMat;
	Avp = new Matrix( "Avp", A, vSolVec, pSolVec, vRHSVec, ops, 1 );

	pDiveMat = new Divergence( "pDive", pressure, velocity, a*dt );
	ops = new Operator*[1];
	ops[0] = pDiveMat;
	Apv = new Matrix( "Apv", A, pSolVec, vSolVec, pRHSVec, ops, 1 );

	pMassMat = new MassMatrix( "pMass", pressure, pressure, 1.0 );
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

void LinearShallowWaterEqn::_Solve( Field* velPrev, Field* presPrev ) {
	Vector*			vSolVec;
	Vector*			pSolVec;
	Vector*			vRHSVec;
	Vector*			pRHSVec;
	FieldRHS*		vRHS;
	FieldRHS*		vRHS2		= NULL;
	GradFieldRHS*		vTopoRHS 	= NULL;
	WindStressRHS*		vWindRHS	= NULL;
	FieldRHS*		pRHS;
	FieldRHS*		pRHS2		= NULL;
	RHSOp**			rhs;
	KSP			ksp, *subksp;
	PC			pc;
	int			nksp;
	IS			uis, vis, pis;
	MatNullSpace		null;
	int			uSize		= velocity->mesh->nVertsTotal - velocity->bcs->size[0];
	int			vSize		= velocity->mesh->nVertsTotal - velocity->bcs->size[1];
	int			pSize 		= pressure->mesh->nVertsTotal - pressure->bcs->size[0];
	double			a		= ( velPrev ) ? 4.0/3.0 : 1.0;
	int			nVelOps 	= 1;
	int			nPresOps	= 1;
	int			velOp_i		= 0;
	int			presOp_i	= 0;
	Field*			pWindStress	= new Field( "pWindStress", pressure->mesh, 1, NULL );

	vSolVec = new Vector( "vSol", velocity, x, NULL, 0 );
	pSolVec = new Vector( "pSol", pressure, x, NULL, 0 );

	if( topography   ) { nVelOps++; }
	if( tau > 1.0e-6 ) { nVelOps++; }
	if( velPrev      ) { nVelOps++; }
	rhs = new RHSOp*[nVelOps];
	vRHS = new FieldRHS( "vRHS", velocity->mesh, a, velocity );
	rhs[velOp_i] = vRHS;
	if( velPrev ) {
		velOp_i++;
		vRHS2 = new FieldRHS( "vRHS2", velocity->mesh, -1.0/3.0, velPrev );
		rhs[velOp_i] = vRHS2;
	}
	if( topography ) {
		velOp_i++;
		vTopoRHS = new GradFieldRHS( "vTopoRHS", velocity->mesh, -dt, topography );
		rhs[velOp_i] = vTopoRHS;
	}
	if( tau > 1.0e-6 ) {
		velOp_i++;
		vWindRHS = new WindStressRHS( "vWindRHS", velocity->mesh, dt*tau, pWindStress, M_PI, 1.0 );
		rhs[velOp_i] = vWindRHS;
	}
	vRHSVec = new Vector( "vRHS", velocity, f, rhs, nVelOps );

	if( presPrev ) { nPresOps++; }
	rhs = new RHSOp*[nPresOps];
	pRHS = new FieldRHS( "pRHS", pressure->mesh, a, pressure );
	rhs[presOp_i] = pRHS;
	if( presPrev ) {
		presOp_i++;
		pRHS2 = new FieldRHS( "pRHS2", pressure->mesh, -1.0/3.0, presPrev );
		rhs[presOp_i] = pRHS2;
	}
	pRHSVec = new Vector( "pRHS", pressure, f, rhs, nPresOps );

	VecZeroEntries( f );
	vRHSVec->Assemble();
	pRHSVec->Assemble();
	VecAXPY( f, 1.0, b );

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
	delete pWindStress;
}

void LinearShallowWaterEqn::CalcEnergetics( int timeStep, double* KE, double* PE, double* visc, double* fric, double* wind ) {
	double detJac, weight, *coord, v[2], p, gCoord[2], **d2v;
	ofstream file;
	char filename[40];

	*KE = *PE = *visc = *fric = *wind = 0.0;

	d2v = new double*[2];
	d2v[0] = new double[2];
	d2v[1] = new double[2];

	for( int el_i = 0; el_i < velocity->mesh->nElsTotal; el_i++ ) {
		for( int pt_i = 0; pt_i < velocity->mesh->el->nPoints; pt_i++ ) {
			coord  = velocity->mesh->el->quadPts[pt_i]->coord;
			weight = velocity->mesh->el->quadPts[pt_i]->weight;
			detJac = velocity->mesh->DetJac( el_i, pt_i );
			velocity->InterpLocal( el_i, coord, v );
			pressure->InterpLocal( el_i, coord, &p );
			velocity->mesh->LocalToGlobal( coord, el_i, gCoord );
			*KE   += detJac*weight*0.5*( 1.0 + p )*( v[0]*v[0] + v[1]*v[1] );
			*PE   += detJac*weight*p*p;
			/* the second derivative operator is a little buggy, leaves an artifact on the element boundaries */
			//*visc += detJac*weight*( 1.0 + p )*( v[0]*( d2v[0][0] + d2v[0][1] ) + v[1]*( d2v[1][0] + d2v[1][1] ) )/detJac;
			*fric += detJac*weight*( 1.0 + p )*( v[0]*v[0] + v[1]*v[1] );
			*wind += detJac*weight*v[0]*cos( M_PI*gCoord[1] );
		}
	}

	delete[] d2v[0];
	delete[] d2v[1];
	delete[] d2v;

	sprintf( filename, "energetics.%.5u.sw", timeStep );
	file.open( filename );
	file << "*****Energetics*****" << endl;
	file << "KE:\t" << *KE << endl;
	file << "PE:\t" << *PE << endl;
	file << "*****Energy Flux****" << endl;
	//file << "viscosity:\t"   << *visc << endl;
	file << "friction:\t"    << *fric << endl;
	file << "wind stress:\t" << *wind << endl;
	file.close();
}
