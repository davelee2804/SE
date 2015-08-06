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

#define EPS 0.02
#define L0 0.0 		// solitary wave number
#define OMEGA 0.01 	//
#define ALPHA +0.35436	// nls eqn dispersive coefficient
#define GAMMA +0.81838	// nls eqn non linear coefficient
#define DELTA -0.13129	// nls eqn first order coefficient
#define KAPPA -0.17484	// nls eqn zeroth order coefficient
#define A0 sqrt( 2*OMEGA/GAMMA )
#define B0 sqrt( OMEGA/ALPHA )
#define C0 0.0
#define NL 12		// no. waves with wavenumber l in domain
#define K0 0.75		// base wave number
#define CP -1.14826389894846	// phase speed
#define F1 0.5		// L^2/Ld^2
#define F2 EPS*F1
#define BETA 2.0	// planetary vorticity gradient
#define CG -0.74503	// group velocity
#define LEN (NL/K0)

double* Ay;
double* dAy;
double* d2Ay;
double* Cy;
double* dCy;
double* d2Cy;

void LoadYStruct( Field* phi1, Field* phi2 ) {
	ifstream	file;

	Ay   = new double[phi1->mesh->nVerts[1]];
	dAy  = new double[phi1->mesh->nVerts[1]];
	d2Ay = new double[phi1->mesh->nVerts[1]];
	Cy   = new double[phi2->mesh->nVerts[1]];
	dCy  = new double[phi2->mesh->nVerts[1]];
	d2Cy = new double[phi2->mesh->nVerts[1]];

	file.open( "Anew.txt" );
	for( int y_i = 0; y_i < phi1->mesh->nVerts[1]; y_i++ ) { 
		file >> Ay[y_i];
	}
	file.close();
	file.open( "dAnew.txt" );
	for( int y_i = 0; y_i < phi1->mesh->nVerts[1]; y_i++ ) { 
		file >> dAy[y_i];
	}
	file.close();
	file.open( "d2Anew.txt" );
	for( int y_i = 0; y_i < phi1->mesh->nVerts[1]; y_i++ ) { 
		file >> d2Ay[y_i];
	}
	file.close();

	file.open( "Cnew.txt" );
	for( int y_i = 0; y_i < phi2->mesh->nVerts[1]; y_i++ ) { 
		file >> Cy[y_i];
	}
	file.close();
	file.open( "dCnew.txt" );
	for( int y_i = 0; y_i < phi2->mesh->nVerts[1]; y_i++ ) { 
		file >> dCy[y_i];
	}
	file.close();
	file.open( "d2Cnew.txt" );
	for( int y_i = 0; y_i < phi2->mesh->nVerts[1]; y_i++ ) { 
		file >> d2Cy[y_i];
	}
	file.close();
}

void GenAnalytic( Field* phi1, Field* phi2, Field* omega1, Field* omega2, double time, bool wShear ) {
	int 	topo[2];
	double	X0, X1, y, theta, zeta;

	for( int vert_i = 0; vert_i < phi1->mesh->nVertsTotal; vert_i++ ) {
		phi1->mesh->IndexToTopo( vert_i, topo );
		X0 = phi1->mesh->verts[vert_i][0];
		X1 = phi1->mesh->verts[vert_i][0] - CG*time;
		y  = phi1->mesh->verts[vert_i][1];
		theta = K0*X0 + L0*X1 - (OMEGA + DELTA*L0 + CP*K0)*time;
		zeta = B0*(X1 - DELTA*time - C0*time);

		phi1->vals[vert_i][0] = EPS*Ay[topo[1]]*2.0*A0*cos( theta )/cosh( zeta );
		if( wShear ) {
			phi1->vals[vert_i][0] -= 0.5*y*y/M_PI;
		}
		phi2->vals[vert_i][0] = EPS*EPS*Cy[topo[1]]*2.0*A0*cos( theta )/cosh( zeta );
		if( wShear ) {
			phi2->vals[vert_i][0] -= EPS*0.5*y*y/M_PI;
		}

		omega1->vals[vert_i][0] = EPS*d2Ay[topo[1]]*2.0*A0*cos( theta )/cosh( zeta ) + 
			EPS*Ay[topo[1]]*2.0*A0*( 2.0*B0*(K0 + L0)*sin( theta )*tanh( zeta ) + 
						 (B0*B0 - (K0 + L0)*(K0 - L0))*cos( theta ) - 
						 B0*B0*cos( theta )/( cosh( zeta )*cosh( zeta ) ) )/cosh( zeta );
		if( wShear ) {
			omega1->vals[vert_i][0] -= 1.0/M_PI;
		}
		omega2->vals[vert_i][0] = EPS*EPS*d2Cy[topo[1]]*2.0*A0*cos( theta )/cosh( zeta ) + 
			EPS*EPS*Cy[topo[1]]*2.0*A0*( 2.0*B0*(K0 + L0)*sin( theta )*tanh( zeta ) + 
						 (B0*B0 - (K0 + L0)*(K0 - L0))*cos( theta ) - 
						 B0*B0*cos( theta )/( cosh( zeta )*cosh( zeta ) ) )/cosh( zeta );
		if( wShear ) {
			omega2->vals[vert_i][0] -= EPS*1.0/M_PI;
		}

		//phi1->vals[vert_i][0] *= cos( KAPPA*time );
		//phi2->vals[vert_i][0] *= cos( KAPPA*time );
		//omega1->vals[vert_i][0] *= cos( KAPPA*time );
		//omega2->vals[vert_i][0] *= cos( KAPPA*time );
	}
}

void GenVelocity( Field* vel1, Field* vel2, double time, bool wShear ) {
	int		topo[2];
	double		X0, X1, y, theta, zeta;

	for( int vert_i = 0; vert_i < vel1->mesh->nVertsTotal; vert_i++ ) {
		vel1->mesh->IndexToTopo( vert_i, topo );
		X0 = vel1->mesh->verts[vert_i][0];
		X1 = vel1->mesh->verts[vert_i][0] - CG*time;
		y  = vel1->mesh->verts[vert_i][1];
		theta = K0*X0 + L0*X1 - (OMEGA + DELTA*L0 + CP*K0)*time;
		zeta = B0*(X1 - DELTA*time - C0*time);

		vel1->vals[vert_i][0] = -EPS*dAy[topo[1]]*2.0*A0*cos( theta )/cosh( zeta );
		if( wShear ) {
			vel1->vals[vert_i][0] += y/M_PI;
		}
		vel1->vals[vert_i][1] = -EPS*2.0*Ay[topo[1]]*A0*( 
					(K0 + L0)*sin( theta ) + B0*cos( theta )*tanh( zeta ) )/cosh( zeta );

		vel2->vals[vert_i][0] = -EPS*EPS*dCy[topo[1]]*2.0*A0*cos( theta )/cosh( zeta );
		if( wShear ) {
			vel2->vals[vert_i][0] += EPS*y/M_PI;
		}
		vel2->vals[vert_i][1] = -EPS*EPS*2.0*Cy[topo[1]]*A0*( 
					(K0 + L0)*sin( theta ) + B0*cos( theta )*tanh( zeta ) )/cosh( zeta );

		//vel1->vals[vert_i][0] *= cos( KAPPA*time );
		//vel1->vals[vert_i][1] *= cos( KAPPA*time );
		//vel2->vals[vert_i][0] *= cos( KAPPA*time );
		//vel2->vals[vert_i][1] *= cos( KAPPA*time );
	}
}

/* batoclinic and barotropic fields */
void CalcTauAndChi( Field* psi1, Field* psi2, Field* tau, Field* chi ) {
	for( int node_i = 0; node_i < tau->mesh->nVertsTotal; node_i++ ) {
		tau->vals[node_i][0] = 0.5*(psi1->vals[node_i][0] - psi2->vals[node_i][0]);
		chi->vals[node_i][0] = 0.5*(psi1->vals[node_i][0] + psi2->vals[node_i][0]);
	}
}

double sf0( double* x ) { return 0.0; }
double sf1( double* x ) { return x[1]/M_PI; }
double sf2( double* x ) { return EPS*x[1]/M_PI; }

double height( double* x ) { return 0.1*cos( 2.5*x[1] )*cos( 6*x[0]/LEN ); }

int main( int argc, char** argv ) {
	char		tag[6]		= "petsc";
	int		nx[2]		= { 192, 12 };
	bool		periodic[2]	= { true, false };
	double		min[2]		= { -LEN*M_PI, -M_PI };
	double		max[2]		= { +LEN*M_PI, +M_PI };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", 5, min, max, periodic );
	BCs*		bcs		= new BCs( true, true, false, false, mesh, 0 );
	Field*		phi1		= new Field( "phi-1-curr", mesh, 1, bcs );
	Field*		phi2		= new Field( "phi-2-curr", mesh, 1, bcs );
	Field*		omega1		= new Field( "omega-1-curr", mesh, 1, NULL );
	Field*		omega2		= new Field( "omega-2-curr", mesh, 1, NULL );
	Field*		phi1Prev	= new Field( "phi-1-prev", mesh, 1, bcs );
	Field*		phi2Prev	= new Field( "phi-2-prev", mesh, 1, bcs );
	Field*		omega1Prev	= new Field( "omega-1-prev", mesh, 1, NULL );
	Field*		omega2Prev	= new Field( "omega-2-prev", mesh, 1, NULL );
	Field*		phi1Anal	= new Field( "phi-1-anal", mesh, 1, NULL );
	Field*		phi2Anal	= new Field( "phi-2-anal", mesh, 1, NULL );
	Field*		omega1Anal	= new Field( "omega-1-anal", mesh, 1, NULL );
	Field*		omega2Anal	= new Field( "omega-2-anal", mesh, 1, NULL );
	Field*		phi1MS		= new Field( "phi-1-ms", mesh, 1, NULL );
	Field*		phi2MS		= new Field( "phi-2-ms", mesh, 1, NULL );
	Field*		omega1MS	= new Field( "omega-1-ms", mesh, 1, NULL );
	Field*		omega2MS	= new Field( "omega-2-ms", mesh, 1, NULL );
	Field* 		vel1Anal 	= new Field( "vel-1-anal", mesh, 2, NULL );
	Field* 		vel2Anal 	= new Field( "vel-2-anal", mesh, 2, NULL );
	Field*		tau		= new Field( "tau", mesh, 1, NULL );
	Field*		chi		= new Field( "chi", mesh, 1, NULL );
	Field*		tauAnal		= new Field( "tau-anal", mesh, 1, NULL );
	Field*		chiAnal		= new Field( "chi-anal", mesh, 1, NULL );
	Field*		topo		= new Field( "topography", mesh, 1, NULL );
	Field*		fields[16];
	double		dt		= 0.01;
	int		timeStep	= atoi( argv[1] );
	double		time		= timeStep*dt;
	int		dumpEvery	= 1;
	int		nTimeSteps	= 400;
	TwoLayerQGEqn*	qg;
	double		y;

	PetscInitialize( &argc, &argv, (char)0, tag );

	LoadYStruct( phi1, phi2 );

	cout << "A0:  " << A0 << endl;
	cout << "Len: " << LEN << endl;
	cout << "c_p: " << CP << endl;
	cout << "c_g: " << CG << endl;
	cout << "C0:  " << C0 << endl;
	cout << "B0:  " << B0 << endl;

	if( timeStep > 0 ) {
		phi1->Read( timeStep );
		phi2->Read( timeStep );
	}
	else {
		//GenAnalytic( phi1, phi2, omega1, omega2, time, true );
		GenAnalytic( phi1, phi2, omega1, omega2, time, false );
	}
	fields[0]  = phi1MS;
	fields[1]  = phi1Anal;
	fields[2]  = phi2MS;
	fields[3]  = phi2Anal;
	fields[4]  = omega1MS;
	fields[5]  = omega1Anal;
	fields[6]  = omega2MS;
	fields[7]  = omega2Anal;
	fields[9]  = vel1Anal;
	fields[11] = vel2Anal;
	fields[12] = tau;
	fields[13] = tauAnal;
	fields[14] = chi;
	fields[15] = chiAnal;

	mesh->Save();

	topo->SetICFunc( 0, height );
	WriteXDMFHeader( 99999 );
	WriteXDMF( &topo, 1, 99999, 0.0, 0.0 );
	WriteXDMFFooter( 99999 );

	qg = new TwoLayerQGEqn( omega1, phi1, F1, omega2, phi2, F2, BETA, dt, sf1, sf2, topo );

	if( timeStep == 0 ) {
		GenAnalytic( phi1Anal, phi2Anal, omega1Anal, omega2Anal, time, false );
		CalcTauAndChi( phi1Anal, phi2Anal, tauAnal, chiAnal );
		fields[8]  = qg->GenVel( phi1, "vel-1", sf0 );
		fields[10] = qg->GenVel( phi2, "vel-2", sf0 );
		WriteXDMFHeader( timeStep );
		WriteXDMF( fields, 16, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
		delete fields[8];
		delete fields[10];
	}

	phi1Prev->Copy( phi1 );
	phi2Prev->Copy( phi2 );
	omega1Prev->Copy( omega1 );
	omega2Prev->Copy( omega2 );

	time += dt;
	timeStep++;

	//GenAnalytic( phi1, phi2, omega1, omega2, time, true );
	GenAnalytic( phi1, phi2, omega1, omega2, time, false );
	CalcTauAndChi( phi1Anal, phi2Anal, tauAnal, chiAnal );
	fields[8]  = qg->GenVel( phi1, "vel-1", sf0 );
	fields[10] = qg->GenVel( phi2, "vel-2", sf0 );
	phi1Prev->Copy( phi1 );
	phi2Prev->Copy( phi2 );
	omega1Prev->Copy( omega1 );
	omega2Prev->Copy( omega2 );
	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 12, timeStep, time, dt );
	WriteXDMFFooter( timeStep );
	delete fields[8];
	delete fields[10];

	for( timeStep++; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		qg->Solve( omega1Prev, phi1Prev, omega2Prev, phi2Prev );
		if( timeStep%dumpEvery == 0 ) {
			GenVelocity( vel1Anal, vel2Anal, time, false );
			fields[8]  = qg->GenVel( phi1, "vel-1", sf0 );
			fields[10] = qg->GenVel( phi2, "vel-2", sf0 );
			for( int i = 0; i < mesh->nVertsTotal; i++ ) {
				y = mesh->verts[i][1];
				//fields[8]->vals[i][0]  -= y/M_PI;
				//fields[10]->vals[i][0] -= EPS*y/M_PI;
				phi1MS->vals[i][0] = phi1->vals[i][0];// + 0.5*y*y/M_PI;
				phi2MS->vals[i][0] = phi2->vals[i][0];// + EPS*0.5*y*y/M_PI;
				omega1MS->vals[i][0] = omega1->vals[i][0];// + 1.0/M_PI;
				omega2MS->vals[i][0] = omega2->vals[i][0];// + EPS*1.0/M_PI;
			}
			
			GenAnalytic( phi1Anal, phi2Anal, omega1Anal, omega2Anal, time, false );
			CalcTauAndChi( phi1, phi2, tau, chi );
			CalcTauAndChi( phi1Anal, phi2Anal, tauAnal, chiAnal );
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 16, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
			delete fields[8];
			delete fields[10];
		}
	}

	delete qg;
	delete mesh;
	delete bcs;
	delete phi1;
	delete phi2;
	delete omega1;
	delete omega2;
	delete phi1Prev;
	delete phi2Prev;
	delete omega1Prev;
	delete omega2Prev;
	delete phi1Anal;
	delete phi2Anal;
	delete omega1Anal;
	delete omega2Anal;
	delete vel1Anal;
	delete vel2Anal;
	delete phi1MS;
	delete phi2MS;
	delete omega1MS;
	delete omega2MS;
	delete tau;
	delete chi;
	delete tauAnal;
	delete chiAnal;
	delete topo;
	delete[] Ay;
	delete[] dAy;
	delete[] d2Ay;
	delete[] Cy;
	delete[] dCy;
	delete[] d2Cy;

	PetscFinalize();

	return EXIT_SUCCESS;
}
