#include <cmath>
#include <fstream>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "LinSys.h"
#include "TwoLayerQGEqn.h"

using namespace std;
using std::string;

void CalcTauAndChi( Field* psi1, Field* psi2, Field* tau, Field* chi ) {
	for( int node_i = 0; node_i < tau->mesh->nVertsTotal; node_i++ ) {
		tau->vals[node_i][0] = 0.5*(psi1->vals[node_i][0] - psi2->vals[node_i][0]);
		chi->vals[node_i][0] = 0.5*(psi1->vals[node_i][0] + psi2->vals[node_i][0]);
	}
}

double* A1;
double* A2;
double* dA1;
double* dA2;
double* d2A1;
double* d2A2;

double* A1_i;
double* A2_i;
double* dA1_i;
double* dA2_i;
double* d2A1_i;
double* d2A2_i;

void LoadYStruct( int ny ) {
	ifstream	file;

	A1   = new double[ny];
	A2   = new double[ny];
	dA1  = new double[ny];
	dA2  = new double[ny];
	d2A1 = new double[ny];
	d2A2 = new double[ny];

	file.open( "A1u.txt" );   for( int y_i = 0; y_i < ny; y_i++ ) { file >> A1[y_i];   } file.close();
	file.open( "A2u.txt" );   for( int y_i = 0; y_i < ny; y_i++ ) { file >> A2[y_i];   } file.close();
	file.open( "dA1u.txt" );  for( int y_i = 0; y_i < ny; y_i++ ) { file >> dA1[y_i];  } file.close();
	file.open( "dA2u.txt" );  for( int y_i = 0; y_i < ny; y_i++ ) { file >> dA2[y_i];  } file.close();
	file.open( "d2A1u.txt" ); for( int y_i = 0; y_i < ny; y_i++ ) { file >> d2A1[y_i]; } file.close();
	file.open( "d2A2u.txt" ); for( int y_i = 0; y_i < ny; y_i++ ) { file >> d2A2[y_i]; } file.close();
}

void LoadYStruct_Imag( int ny ) {
	ifstream	file;

	A1_i   = new double[ny];
	A2_i   = new double[ny];
	dA1_i  = new double[ny];
	dA2_i  = new double[ny];
	d2A1_i = new double[ny];
	d2A2_i = new double[ny];

	file.open( "A1u_i.txt" );   for( int y_i = 0; y_i < ny; y_i++ ) { file >> A1_i[y_i];   } file.close();
	file.open( "A2u_i.txt" );   for( int y_i = 0; y_i < ny; y_i++ ) { file >> A2_i[y_i];   } file.close();
	file.open( "dA1u_i.txt" );  for( int y_i = 0; y_i < ny; y_i++ ) { file >> dA1_i[y_i];  } file.close();
	file.open( "dA2u_i.txt" );  for( int y_i = 0; y_i < ny; y_i++ ) { file >> dA2_i[y_i];  } file.close();
	file.open( "d2A1u_i.txt" ); for( int y_i = 0; y_i < ny; y_i++ ) { file >> d2A1_i[y_i]; } file.close();
	file.open( "d2A2u_i.txt" ); for( int y_i = 0; y_i < ny; y_i++ ) { file >> d2A2_i[y_i]; } file.close();
}
//double uTopIC( double* x ) { return exp( -0.1*(x[1]-10.5)*(x[1]-10.5) )*cos( K0*x[0] ); }
//double omegaIC( double* x ) { return 0.2*(x[1] - 10.5)*exp( -0.1*(x[1]-10.5)*(x[1]-10.5) ); }
double uTopIC( double* x )  { return 1.0/cosh( 0.5*x[1] )/cosh( 0.5*x[1] ); }
double psi1IC( double* x )  { return -2.0*tanh( 0.5*x[1] ); }
double omegaIC( double* x ) { return tanh( 0.5*x[1] )/( cosh( 0.5*x[1] )*cosh( 0.5*x[1] ) ); }

void ShearComponents( Field* phi, Field* omega, Field* vel, bool add ) {
	double	fac	= ( add ) ? +1.0 : -1.0;
	double*	x;

	for( int i = 0; i < phi->mesh->nVertsTotal; i++ ) {
		x = phi->mesh->verts[i];

		//phi->vals[i][0]   += fac*mean->vals[i][0];
		//omega->vals[i][0] += fac*0.2*(x[1]-10.5)*exp( -0.1*(x[1]-10.5)*(x[1]-10.5) );
		//vel->vals[i][0]   += fac*exp( -0.1*(x[1]-10.5)*(x[1]-10.5) );
		phi->vals[i][0]   -= fac*2.0*tanh( 0.5*x[1] );
		omega->vals[i][0] += fac*tanh( 0.5*x[1] )/( cosh( 0.5*x[1] )*cosh( 0.5*x[1] ) );
		if( vel ) {
			vel->vals[i][0]   += fac/( cosh( 0.5*x[1] )*cosh( 0.5*x[1] ) );
		}
	}
}

#define EPS 0.02
#define K0 0.8
#define CPR  +0.157358 	// unstable
#define CPI  +0.106444	// unstable
//#define CPR  -0.746001 	// stable
//#define CPI  +0.0	// stable
#define BETA 0.2
//#define BETA 0.5

/*
void FirstOrdComponents( Field* phi1, Field* omega1, Field* vel1, Field* phi2, Field* omega2, Field* vel2, double time ) {
	double*	x;
	int	topo[2], yi;
	double	amp	= exp( K0*CPI*time );

	for( int i = 0; i < phi1->mesh->nVertsTotal; i++ ) {
		x = phi1->mesh->verts[i];
		phi1->mesh->IndexToTopo( i, topo );
		yi = topo[1];

		phi1->vals[i][0]   = +amp*EPS*A1[yi]*cos( K0*x[0] - CPR*K0*time );
		omega1->vals[i][0] = +amp*EPS*(d2A1[yi] - K0*K0*A1[yi])*cos( K0*x[0] - CPR*K0*time );
		if( vel1 ) {
			vel1->vals[i][0] = -amp*EPS*dA1[yi]*cos( K0*x[0] - CPR*K0*time );
			vel1->vals[i][1] = -amp*EPS*K0*A1[yi]*sin( K0*x[0] - CPR*K0*time );
		}

		phi2->vals[i][0]   = +amp*EPS*A2[yi]*cos( K0*x[0] - CPR*K0*time );
		omega2->vals[i][0] = +amp*EPS*(d2A2[yi] - K0*K0*A2[yi])*cos( K0*x[0] - CPR*K0*time );
		if( vel2 ) {
			vel2->vals[i][0] = -amp*EPS*dA2[yi]*cos( K0*x[0] - CPR*K0*time );
			vel2->vals[i][1] = -amp*EPS*K0*A2[yi]*sin( K0*x[0] - CPR*K0*time );
		}
	}
}
*/

void FirstOrdComponents( Field* phi1, Field* omega1, Field* vel1, Field* phi2, Field* omega2, Field* vel2, double time ) {
	double*	x;
	int	topo[2], yi;
	double	amp	= exp( K0*CPI*time );

	for( int i = 0; i < phi1->mesh->nVertsTotal; i++ ) {
		x = phi1->mesh->verts[i];
		phi1->mesh->IndexToTopo( i, topo );
		yi = topo[1];

		phi1->vals[i][0]   = +amp*EPS*A1[yi]*cos( K0*x[0] - CPR*K0*time ) - amp*EPS*A1_i[yi]*sin( K0*x[0] - CPR*K0*time );
		omega1->vals[i][0] = +amp*EPS*(d2A1[yi] - K0*K0*A1[yi])*cos( K0*x[0] - CPR*K0*time ) - amp*EPS*(d2A1_i[yi] - K0*K0*A1_i[yi])*sin( K0*x[0] - CPR*K0*time );
		if( vel1 ) {
			vel1->vals[i][0] = -amp*EPS*dA1[yi]*cos( K0*x[0] - CPR*K0*time ) + amp*EPS*dA1_i[yi]*sin( K0*x[0] - CPR*K0*time );
			vel1->vals[i][1] = -amp*EPS*K0*A1[yi]*sin( K0*x[0] - CPR*K0*time ) - amp*EPS*K0*A1_i[yi]*cos( K0*x[0] - CPR*K0*time );
		}

		phi2->vals[i][0]   = +amp*EPS*A2[yi]*cos( K0*x[0] - CPR*K0*time ) - amp*EPS*A2_i[yi]*sin( K0*x[0] - CPR*K0*time );
		omega2->vals[i][0] = +amp*EPS*(d2A2[yi] - K0*K0*A2[yi])*cos( K0*x[0] - CPR*K0*time ) - amp*EPS*(d2A2_i[yi] - K0*K0*A2_i[yi])*sin( K0*x[0] - CPR*K0*time );
		if( vel2 ) {
			vel2->vals[i][0] = -amp*EPS*dA2[yi]*cos( K0*x[0] - CPR*K0*time ) + amp*EPS*dA2_i[yi]*sin( K0*x[0] - CPR*K0*time );
			vel2->vals[i][1] = -amp*EPS*K0*A2[yi]*sin( K0*x[0] - CPR*K0*time ) - amp*EPS*K0*A2_i[yi]*cos( K0*x[0] - CPR*K0*time );
		}
	}
}

#define NK 8
#define LEN (2.0*M_PI*NK/K0)

void TestPsi( Field* psi ) {
	Field*		fields[2];
	Field*		uxa	= new Field( "uxa", psi->mesh, 1, NULL );
	Field*		uxn	= new Field( "uxn", psi->mesh, 1, NULL );
	double*		x;
	double**	du	= new double*[1];

	du[0] = new double[2];
	fields[0] = uxa;
	fields[1] = uxn;

	for( int i = 0; i < psi->mesh->nVertsTotal; i++ ) {
		x = psi->mesh->verts[i];
		psi->InterpDerivsGlobal( x, du );
		uxa->vals[i][0] = 1.0/cosh(0.5*x[1])/cosh(0.5*x[1]);
		uxn->vals[i][0] = -du[0][1];
	}
	WriteXDMFHeader( 99999 );
	WriteXDMF( fields, 2, 99999, 0.0, 0.0 );
	WriteXDMFFooter( 99999 );

	delete[] du[0];
	delete[] du;
	delete uxa;
	delete uxn;
}

int main( int argc, char** argv ) {
	char		tag[6]		= "petsc";
	int		nx[2]		= { 100, 32 };
	bool		periodic[2]	= { true, false };
	double		min[2]		= { 0.0, -20.0 };
	double		max[2]		= { LEN, +20.0 };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", 4, min, max, periodic );
	BCs*		bcs		= new BCs( true, true, false, false, mesh, 0 );
	Field*		psi1		= new Field( "psi-1-curr", mesh, 1, bcs );
	Field*		psi2		= new Field( "psi-2-curr", mesh, 1, bcs );
	Field*		omega1		= new Field( "omega-1-curr", mesh, 1, bcs );
	Field*		omega2		= new Field( "omega-2-curr", mesh, 1, bcs );
	Field*		psi1Prev	= new Field( "psi-1-prev", mesh, 1, NULL );
	Field*		psi2Prev	= new Field( "psi-2-prev", mesh, 1, NULL );
	Field*		omega1Prev	= new Field( "omega-1-prev", mesh, 1, NULL );
	Field*		omega2Prev	= new Field( "omega-2-prev", mesh, 1, NULL );
	Field*		tau		= new Field( "tau", mesh, 1, NULL );
	Field*		chi		= new Field( "chi", mesh, 1, NULL );
	Field*		tauAnal		= new Field( "tau-anal", mesh, 1, NULL );
	Field*		chiAnal		= new Field( "chi-anal", mesh, 1, NULL );
	Field*		tau_e		= new Field( "tau-e", mesh, 1, bcs );
	Field*		psi1Anal	= new Field( "psi-1-anal", mesh, 1, NULL );
	Field*		psi2Anal	= new Field( "psi-2-anal", mesh, 1, NULL );
	Field*		omega1Anal	= new Field( "omega-1-anal", mesh, 1, NULL );
	Field*		omega2Anal	= new Field( "omega-2-anal", mesh, 1, NULL );
	Field*		vel1Anal	= new Field( "vel-1-anal", mesh, 2, NULL );
	Field*		vel2Anal	= new Field( "vel-2-anal", mesh, 2, NULL );
	Field*		fields[16];
	double		dt		= 0.01;
	int		timeStep	= atoi( argv[1] );
	double		time		= timeStep*dt;
	int		dumpEvery	= 1;
	int		nTimeSteps	= 400;
	TwoLayerQGEqn*	qg;

	PetscInitialize( &argc, &argv, (char)0, tag );

	LoadYStruct( mesh->nVerts[1] );
	LoadYStruct_Imag( mesh->nVerts[1] );

	if( timeStep > 0 ) {
		psi1->Read( timeStep );
		psi2->Read( timeStep );
		omega1->Read( timeStep );
		omega2->Read( timeStep );
	}
	else {
		omega1->SetICFunc( 0, omegaIC );
		psi1->SetICFunc( 0, psi1IC );
		tau_e->SetICFunc( 0, psi1IC );
		FirstOrdComponents( psi1Anal, omega1Anal, NULL, psi2Anal, omega2Anal, NULL, 0.0 );
		for( int i = 0; i < mesh->nVertsTotal; i++ ) {
			psi1->vals[i][0]   += psi1Anal->vals[i][0];
			psi2->vals[i][0]   += psi2Anal->vals[i][0];
			omega1->vals[i][0] += omega1Anal->vals[i][0];
			omega2->vals[i][0] += omega2Anal->vals[i][0];
		}
	}
	fields[0]  = psi1;
	fields[1]  = psi1Anal;
	fields[2]  = psi2;
	fields[3]  = psi2Anal;
	fields[4]  = omega1;
	fields[5]  = omega1Anal;
	fields[6]  = omega2;
	fields[7]  = omega2Anal;
	fields[9]  = vel1Anal;
	fields[11] = vel2Anal;
	fields[12] = chi;
	fields[13] = chiAnal;
	fields[14] = tau;
	fields[15] = tauAnal;

	mesh->Save();

	qg = new TwoLayerQGEqn( psi1, omega1, psi2, omega2, 0.5, 0.5, BETA, 0.0, 0.0, 0.0001, tau_e, dt );

	if( timeStep == 0 ) {
		fields[8]  = qg->GenVel( psi1, "vel-1" );
		fields[10] = qg->GenVel( psi2, "vel-2" );
		ShearComponents( psi1, omega1, fields[8], false );
		CalcTauAndChi( psi1, psi2, tau, chi );
		CalcTauAndChi( psi1Anal, psi2Anal, tauAnal, chiAnal );
		FirstOrdComponents( psi1Anal, omega1Anal, vel1Anal, psi2Anal, omega2Anal, vel2Anal, 0.0 );
		WriteXDMFHeader( timeStep );
		WriteXDMF( fields, 16, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
		ShearComponents( psi1, omega1, fields[8], true );
		delete fields[8];
		delete fields[10];
	}
	psi1Prev->Copy( psi1 );
	psi2Prev->Copy( psi2 );
	omega1Prev->Copy( omega1 );
	omega2Prev->Copy( omega2 );

	time += dt;
	timeStep++;
	qg->Solve( NULL, NULL, NULL, NULL );
	/*{
		omega1->SetICFunc( 0, omegaIC );
		psi1->SetICFunc( 0, psi1IC );
		FirstOrdComponents( psi1Anal, omega1Anal, NULL, psi2Anal, omega2Anal, NULL, time );
		for( int i = 0; i < mesh->nVertsTotal; i++ ) {
			psi1->vals[i][0]   += psi1Anal->vals[i][0];
			psi2->vals[i][0]   += psi2Anal->vals[i][0];
			omega1->vals[i][0] += omega1Anal->vals[i][0];
			omega2->vals[i][0] += omega2Anal->vals[i][0];
		}
	}*/

	fields[8]  = qg->GenVel( psi1, "vel-1" );
	fields[10] = qg->GenVel( psi2, "vel-2" );
	ShearComponents( psi1, omega1, fields[8], false );
	CalcTauAndChi( psi1, psi2, tau, chi );
	CalcTauAndChi( psi1Anal, psi2Anal, tauAnal, chiAnal );
	FirstOrdComponents( psi1Anal, omega1Anal, vel1Anal, psi2Anal, omega2Anal, vel2Anal, time );
	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 16, timeStep, time, dt );
	WriteXDMFFooter( timeStep );
	ShearComponents( psi1, omega1, fields[8], true );
	delete fields[8];
	delete fields[10];

	for( timeStep++; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		qg->Solve( psi1Prev, omega1Prev, psi2Prev, omega2Prev );
		if( timeStep%dumpEvery == 0 ) {
			fields[8]  = qg->GenVel( psi1, "vel-1" );
			fields[10] = qg->GenVel( psi2, "vel-2" );
			ShearComponents( psi1, omega1, fields[8], false );
			CalcTauAndChi( psi1, psi2, tau, chi );
			CalcTauAndChi( psi1Anal, psi2Anal, tauAnal, chiAnal );
			FirstOrdComponents( psi1Anal, omega1Anal, vel1Anal, psi2Anal, omega2Anal, vel2Anal, time );
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 16, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
			ShearComponents( psi1, omega1, fields[8], true );
			delete fields[8];
			delete fields[10];
		}
	}

	delete qg;
	delete mesh;
	delete bcs;
	delete psi1;
	delete psi2;
	delete omega1;
	delete omega2;
	delete psi1Prev;
	delete psi2Prev;
	delete omega1Prev;
	delete omega2Prev;
	delete tau;
	delete chi;
	delete tauAnal;
	delete chiAnal;
	delete tau_e;
	delete psi1Anal;
	delete psi2Anal;
	delete omega1Anal;
	delete omega2Anal;
	delete vel1Anal;
	delete vel2Anal;
	delete[] A1;
	delete[] A2;
	delete[] dA1;
	delete[] dA2;
	delete[] d2A1;
	delete[] d2A2;

	PetscFinalize();

	return EXIT_SUCCESS;
}
