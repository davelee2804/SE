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
#include "ShallowWaterEqnSC.h"

using namespace std;
using std::string;

ShallowWaterEqnSC::ShallowWaterEqnSC( Field* _vel, Field* _eta, double _dt, ShallowWaterParams* _params, Field* _topo, Field* _velBar, int svvCutoff ) {
	vel 	= _vel;
	eta 	= _eta;
	dt  	= _dt;
	params 	= _params;

	hh = new HelmholtzEqn( vel, vel, 1.0, params->nu, dt, svvCutoff );
	/* assume that the bcs are homogenous, otherwise they're being applied twice!! */
	hh->Assemble();

	sc = NULL;
}

ShallowWaterEqnSC::~ShallowWaterEqnSC() {
	MatDestroy( Kinv );
	MatDestroy( G );
	MatDestroy( D );
	MatDestroy( C );
	VecDestroy( b );
	delete hh;
	delete sc;
}

void ShallowWaterEqnSC::Solve( int order, Field* velRhsField ) {
	Vec		v;
	Vec		e;
	Vec		f;
	Vec		h;
	Vector*		vSol;
	Vector*		eSol;
	Vector*		vRhs;
	Vector*		eRhs;
	DivVelBarRHS*	eDivVelBarRHS	= NULL;
	double		a		= ( order == 2 ) ? dt*2.0/3.0 : dt*1.0;

	CreateVector( vel, &v );
	CreateVector( eta, &e );
	CreateVector( vel, &f );
	CreateVector( eta, &h );

	vSol = new Vector( "v", vel, v, NULL, 0 );
	eSol = new Vector( "e", eta, e, NULL, 0 );

	if( velBar ) { eDivVelBarRHS = new DivVelBarRHS( "eDivVelBarRHS", eta->mesh, a, params->H_i, topo, velBar ); }

	sc->Solve( v, e, f, h );
	
	VecDestroy( v );
	VecDestroy( e );
	VecDestroy( f );
	VecDestroy( h );
	delete vSol;
	delete eSol;
	delete vRhs;
	delete eRhs;
}

void ShallowWaterEqnSC::Assemble( int order ) {
	Matrix* 		vvMat;
	Matrix*			veMat;
	Matrix*			evMat;
	Matrix* 		eeMat;
	BetaInvMatrix*		vBetaInvMat;
	Gradient*		vGradientMat;
	BarotropicDivergence*	eDivergenceMat;
	MassMatrix*		eMassMat;
	DivVelBarMatrix*	eDivVelBarMat	= NULL;
	Operator**		ops;
	Vector*			vSol;
	Vector*			eSol;
	Vector*			vRhs;
	Vec			v;
	Vec			e;
	double			a		= ( order == 2 ) ? 2.0*dt/3.0 : 1.0*dt;

	if( sc ) {
		MatDestroy( Kinv );
		MatDestroy( G );
		MatDestroy( D );
		MatDestroy( C );
		VecDestroy( b );
		delete sc;
	}

	CreateMatrix( vel, vel, &Kinv );
	CreateMatrix( vel, eta, &G );
	CreateMatrix( eta, vel, &D );
	CreateMatrix( eta, eta, &C );
	CreateVector( vel, &b );
	CreateVector( vel, &v );
	CreateVector( eta, &e );

	vSol = new Vector( "v", vel, v, NULL, 0 );
	eSol = new Vector( "e", eta, e, NULL, 0 );
	vRhs = new Vector( "b", vel, b, NULL, 0 );

	vBetaInvMat = new BetaInvMatrix( "vBetaMat", vel, vel, a, params->gamma, params->f0, params->beta );
	ops = new Operator*[1];
	ops[0] = vBetaInvMat;
	vvMat = new Matrix( "Kinv", Kinv, vSol, vSol, vRhs, ops, 1 );
	vvMat->Assemble();

	vGradientMat = new Gradient( "vGradientMat", vel, eta, a*params->g );
	ops = new Operator*[1];
	ops[0] = vGradientMat;
	veMat = new Matrix( "G", G, vSol, eSol, vRhs, ops, 1 );
	veMat->Assemble();

	eDivergenceMat = new BarotropicDivergence( "eDivergenceMat", eta, vel, params->H_i, topo );
	ops = new Operator*[1];
	ops[0] = eDivergenceMat;
	evMat = new Matrix( "D", D, eSol, vSol, NULL, ops, 1 );
	evMat->Assemble();

	eMassMat = new MassMatrix( "eMassMat", eta, eta, 1.0 );
	if( velBar ) { eDivVelBarMat = new DivVelBarMatrix( "eDivVelBarMat", eta, eta, a, velBar ); }
	ops = new Operator*[velBar?2:1];
	ops[0] = eMassMat;
	if( velBar ) { ops[1] = eDivVelBarMat; }
	eeMat = new Matrix( "C", C, eSol, eSol, NULL, ops, velBar?2:1 );
	eeMat->Assemble();

	sc = new SchurComplement( Kinv, G, D, C );

	VecDestroy( v );
	VecDestroy( e );
	delete vSol;
	delete eSol;
	delete vRhs;
	delete vvMat;
	delete veMat;
	delete evMat;
	delete eeMat;
}
