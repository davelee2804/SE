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
#include "QuasiGeostrophicEqn.h"

using namespace std;
using std::string;

#define EPS 0.02
#define L0 0.0 		// solitary wave number
#define OMEGA 0.01 	//
#define ALPHA 0.35436	// nls eqn dispersive coefficient
#define GAMMA 0.81838	// nls eqn non linear coefficient
#define A0 sqrt( 2*OMEGA/GAMMA )
#define B0 sqrt( OMEGA/ALPHA )
#define C0 0.0
#define NL 12		// no. waves with wavenumber l in domain
#define K0 0.75		// base wave number
#define CP -1.14826389894846	// phase speed
#define F 0.5		// L^2/Ld^2
#define BETA 2.0	// planetary vorticity gradient
#define CG -0.74503	// group velocity
#define LEN (NL/K0)

//#define EPS 0.02
//#define L0 0.0 		// solitary wave number
//#define OMEGA 0.01 	//
//#define ALPHA 0.29334	// nls eqn dispersive coefficient
//#define GAMMA 0.20225	// nls eqn non linear coefficient
//#define A0 sqrt( 2*OMEGA/GAMMA )
//#define B0 sqrt( OMEGA/ALPHA )
//#define C0 0.0
//#define NL 12		// no. waves with wavenumber l in domain
//#define K0 0.9		// base wave number
//#define CP -1.289879	// phase speed
//#define F 0.0		// L^2/Ld^2
//#define BETA 2.0	// planetary vorticity gradient
//#define CG -0.45633	// group velocity
//#define LEN (NL/K0)

double* Ay;
double* dAy;
double* d2Ay;

void LoadYStruct( Field* phi ) {
	ifstream	file;

	Ay   = new double[phi->mesh->nVerts[1]];
	dAy  = new double[phi->mesh->nVerts[1]];
	d2Ay = new double[phi->mesh->nVerts[1]];

	file.open( "Anew.txt" );
	for( int y_i = 0; y_i < phi->mesh->nVerts[1]; y_i++ ) { 
		file >> Ay[y_i];
	}
	file.close();
	file.open( "dAnew.txt" );
	for( int y_i = 0; y_i < phi->mesh->nVerts[1]; y_i++ ) { 
		file >> dAy[y_i];
	}
	file.close();
	file.open( "d2Anew.txt" );
	for( int y_i = 0; y_i < phi->mesh->nVerts[1]; y_i++ ) { 
		file >> d2Ay[y_i];
	}
	file.close();
}

void GenAnalytic( Field* phi, Field* omega, double time, bool wShear ) {
	int 	topo[2];
	double	X0, X1, y, theta, zeta;

	for( int vert_i = 0; vert_i < phi->mesh->nVertsTotal; vert_i++ ) {
		phi->mesh->IndexToTopo( vert_i, topo );
		X0 = phi->mesh->verts[vert_i][0];
		X1 = phi->mesh->verts[vert_i][0] - CG*time;
		y  = phi->mesh->verts[vert_i][1];
		theta = K0*X0 + L0*X1 - (CP*K0 - OMEGA)*time;
		zeta = B0*(X1 - C0*time);

		phi->vals[vert_i][0] = EPS*Ay[topo[1]]*2.0*A0*cos( theta )/cosh( zeta );
		if( wShear ) {
			phi->vals[vert_i][0] -= 0.5*y*y/M_PI;
		}

		omega->vals[vert_i][0] = EPS*d2Ay[topo[1]]*2.0*A0*cos( theta )/cosh( zeta ) + 
			EPS*Ay[topo[1]]*2.0*A0*( 2.0*B0*(K0 + L0)*sin( theta )*tanh( zeta ) + 
						 (B0*B0 - (K0 + L0)*(K0 - L0))*cos( theta ) - 
						 B0*B0*cos( theta )/( cosh( zeta )*cosh( zeta ) ) )/cosh( zeta );
		if( wShear ) {
			omega->vals[vert_i][0] -= 1.0/M_PI;
		}
	}
}

void GenVelocity( Field* vela, double time, bool wShear ) {
	int		topo[2];
	double		X0, X1, y, theta, zeta;

	for( int vert_i = 0; vert_i < vela->mesh->nVertsTotal; vert_i++ ) {
		vela->mesh->IndexToTopo( vert_i, topo );
		X0 = vela->mesh->verts[vert_i][0];
		X1 = vela->mesh->verts[vert_i][0] - CG*time;
		y  = vela->mesh->verts[vert_i][1];
		theta = K0*X0 + L0*X1 - (CP*K0 + OMEGA)*time;
		zeta = B0*(X1 - C0*time);

		vela->vals[vert_i][0] = -EPS*dAy[topo[1]]*2.0*A0*cos( theta )/cosh( zeta );
		if( wShear ) {
			vela->vals[vert_i][0] += y/M_PI;
		}
		vela->vals[vert_i][1] = -EPS*2.0*Ay[topo[1]]*A0*( 
					(K0 + L0)*sin( theta ) + B0*cos( theta )*tanh( zeta ) )/cosh( zeta );
	}
}

void WriteAbcissa( Mesh* mesh ) {
	ofstream file;

	file.open( "abcissa.txt" );
	for( int i = 0; i <= mesh->N; i++ ) {
		file << mesh->el->abcissa[i] << endl;
	}
	file.close();
}

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 192, 12 };
	bool			periodic[2]	= { true, false };
	double			min[2]		= { -LEN*M_PI, -M_PI };
	double			max[2]		= { +LEN*M_PI, +M_PI };
	Mesh*			mesh		= new Mesh( "mesh", nx, "legendre", 5, min, max, periodic );
	BCs*			bcs		= new BCs( true, true, false, false, mesh, 0 );
	Field*			phi		= new Field( "phi", mesh, 1, bcs );
	Field*			omega		= new Field( "omega", mesh, 1, NULL );
	Field*			phiPrev		= new Field( "phi-prev", mesh, 1, bcs );
	Field*			omegaPrev	= new Field( "omega-prev", mesh, 1, NULL );
	Field*			phia		= new Field( "phi-anal", mesh, 1, NULL );
	Field*			omegaa		= new Field( "omega-anal", mesh, 1, NULL );
	Field* 			vela 		= new Field( "vel-a", mesh, 2, NULL );
	Field*			phims		= new Field( "phi-ms", mesh, 1, NULL );
	Field*			omegams		= new Field( "omega-ms", mesh, 1, NULL );
	Field* 			veln;
	Field*			fields[8];
	double			dt		= 0.01;
	int			timeStep	= atoi( argv[1] );
	double			time		= timeStep*dt;
	int			dumpEvery	= 1;
	int			nTimeSteps	= 400;
	QuasiGeostrophicEqn*	qg;
	double			y;

	PetscInitialize( &argc, &argv, (char)0, tag );

	LoadYStruct( phi );
	WriteAbcissa( mesh );

	cout << "A0:  " << A0 << endl;
	cout << "Len: " << LEN << endl;
	cout << "c_p: " << CP << endl;
	cout << "c_g: " << CG << endl;
	cout << "C0:  " << C0 << endl;
	cout << "B0:  " << B0 << endl;

	if( timeStep > 0 ) {
		phi->Read( timeStep );
	}
	else {
		GenAnalytic( phi, omega, time, true );
	}
	fields[0] = phi;
	fields[1] = phia;
	fields[2] = omega;
	fields[3] = omegaa;

	//phi->SetBCConst( "bottom", 0, 0.0 );
	//phi->SetBCConst( "top"   , 0, 0.0 );

	mesh->Save();

	if( timeStep == 0 ) {
		GenAnalytic( phia, omegaa, 0.0, false );
		WriteXDMFHeader( timeStep );
		WriteXDMF( fields, 4, timeStep, time, dt );
		WriteXDMFFooter( timeStep );
	}

	phiPrev->Copy( phi );
	omegaPrev->Copy( omega );

	qg = new QuasiGeostrophicEqn( omega, phi, F, BETA, dt );
	//qg->shearFac = 1.0/M_PI;
	qg->shearFac = 0.0;

	time += dt;
	timeStep++;

	GenAnalytic( phi, omega, time, true );
	phia->Copy( phi );
	omegaa->Copy( omega );
	WriteXDMFHeader( timeStep );
	WriteXDMF( fields, 4, timeStep, time, dt );
	WriteXDMFFooter( timeStep );

	for( timeStep++; timeStep <= nTimeSteps; timeStep++ ) {
		time += dt;
		cout << "solving for time step: " << timeStep << "\tdt: " << dt << "\ttime: " << time << endl;
		qg->Solve( omegaPrev, phiPrev );
		if( timeStep%dumpEvery == 0 ) {
			GenVelocity( vela, time, false );
			veln = qg->GenVel( phi );
			for( int i = 0; i < mesh->nVertsTotal; i++ ) {
				y = mesh->verts[i][1];
				veln->vals[i][0] -= y/M_PI;
				phims->vals[i][0] = phi->vals[i][0] + 0.5*y*y/M_PI;
				omegams->vals[i][0] = omega->vals[i][0] + 1.0/M_PI;
			}
			fields[4] = vela;
			fields[5] = veln;
			fields[6] = phims;
			fields[7] = omegams;
			
			GenAnalytic( phia, omegaa, time, false );
			WriteXDMFHeader( timeStep );
			WriteXDMF( fields, 8, timeStep, time, dt );
			WriteXDMFFooter( timeStep );
			WriteXDMFTemporal( timeStep, dumpEvery );
		}
	}

	delete qg;
	delete mesh;
	delete bcs;
	delete phi;
	delete omega;
	delete phiPrev;
	delete omegaPrev;
	delete phia;
	delete omegaa;
	delete phims;
	delete omegams;
	delete[] Ay;
	delete[] dAy;
	delete[] d2Ay;

	PetscFinalize();

	return EXIT_SUCCESS;
}
