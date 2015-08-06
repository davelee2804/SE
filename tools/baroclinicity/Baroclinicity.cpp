#include <cmath>
#include <iostream>
#include <fstream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"

using namespace std;
using std::string;

int main( int argc, char** argv ) {
	char			tag[6]		= "petsc";
	int			nx[2]		= { 18, 18 };
	bool			periodic[2]	= { false, false };
	double			min[2]		= { 0.0, 0.0 };
	double			max[2]		= { 1.0, 1.0 };
	Mesh*			vMesh		= new Mesh( "v-mesh", nx, "legendre", 5, min, max, periodic );
	Mesh*			hMesh		= new Mesh( "h-mesh", nx, "legendre", 3, min, max, periodic );
	Field*			pres		= new Field( "pres-surf", vMesh, 1, NULL );
	Field*			eta		= new Field( "eta-int", hMesh, 1, NULL );
	Field*			baro		= new Field( "pres-int", hMesh, 1, NULL );
	Field*			velTop		= new Field( "vel-top", vMesh, 2, NULL );
	Field*			velBot		= new Field( "vel-bot", vMesh, 2, NULL );
	Field*			velBar		= new Field( "vel-bar", vMesh, 2, NULL );
	int			skip		= atoi( argv[1] );
	int			end		= atoi( argv[2] );
	double			dt		= 0.2/(3.0*nx[0]);
	double			pFac		= 1.0;
	double			hFac		= 20.0;
	Field*			vFields[4];
	Field*			hFields[2];
	double			p, e;

	PetscInitialize( &argc, &argv, (char)0, tag );

	vFields[0] = velTop;
	vFields[1] = velBot;
	vFields[2] = pres;
	vFields[3] = velBar;
	hFields[0] = eta;
	hFields[1] = baro;

	for( int step_i = skip; step_i <= end; step_i += skip ) {
		pres->Read( step_i );
		eta->Read( step_i );
		velTop->Read( step_i );
		velBot->Read( step_i );
		velBar->Read( step_i );

		for( int node_i = 0; node_i < hMesh->nVertsTotal; node_i++ ) {
			pres->InterpGlobal( hMesh->verts[node_i], &p );
			e = eta->vals[node_i][0];
			baro->vals[node_i][0] = pFac*p + hFac*e;
		}

		WriteXDMFHeader( step_i );
		WriteXDMF( hFields, 2, step_i, step_i*dt, dt );
		WriteXDMF( vFields, 4, step_i, step_i*dt, dt );
		WriteXDMFFooter( step_i );
	}
	WriteXDMFTemporal( end, skip );

	delete vMesh;
	delete hMesh;
	delete pres;
	delete eta;
	delete baro;
	delete velTop;
	delete velBot;
	delete velBar;

	PetscFinalize();

	return EXIT_SUCCESS;
}
