#include <iostream>
#include <string>

#include <petsc.h>

#include "Base.h"
#include "LinSys.h"

using namespace std;
using std::string;

#define ELORD 16

int main( int argc, char** argv ) {
	char		tag[]		= "petsc";
	double 		min[2]		= { -10.0, -10.0 };
	double 		max[2]		= { +10.0, +10.0 };
	int		nx[2]		= { 4, 4 };
	bool		periodic[2]	= { false, false };
	Mesh*		mesh		= new Mesh( "mesh", nx, "legendre", ELORD, min, max, periodic );
	BCs*		none		= new BCs( true, true, true, true, mesh, 0 );
	Field*		phiClean	= new Field( "phi-clean", mesh, 1, none );
	Field*		phiModal	= new Field( "phi-modal", mesh, 1, none );
	Field*		phiNodal	= new Field( "phi-nodal", mesh, 1, none );
	Vec		mVec;
	Vec		nVec;
	Vector*		phiModalVec;
	Vector*		phiNodalVec;
	int*		nodes;
	Mat		M2N;
	Field*		fields[2];
	HelmholtzEqn*	he;

	PetscInitialize( &argc, &argv, (char)0, tag );

	/* set up the initial field with the high frequency component */
	for( int el_i = 0; el_i < mesh->nElsTotal; el_i++ ) {
		nodes = mesh->ElNodes( el_i );
		phiModal->vals[nodes[2*(ELORD+1)+2]][0] = 16.0;
		phiModal->vals[nodes[(ELORD-2)*(ELORD+1)+(ELORD-2)]][0] = 16.0;
	}
	CreateModalToNodalMatrix( phiModal, &M2N );
	MatAssemblyBegin( M2N, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( M2N, MAT_FINAL_ASSEMBLY );

	CreateVector( phiModal, &mVec );
	CreateVector( phiNodal, &nVec );
	
	phiModalVec = new Vector( "phi-modal-vec", phiModal, mVec, NULL, 0 );
	phiNodalVec = new Vector( "phi-nodal-vec", phiNodal, nVec, NULL, 0 );

	phiModalVec->UpdateVec();
	MatMult( M2N, mVec, nVec );
	phiNodalVec->UpdateField();

	/* apply the filter */
	//Filter( phiNodal, phiClean );
	he = new HelmholtzEqn( phiClean, phiNodal, 1.0, 0.0, 1.0, ELORD/2 );
	he->Solve( "helm_" );

	fields[0] = phiNodal;
	fields[1] = phiClean;
	mesh->Save();
	WriteXDMFHeader( 0 );
	WriteXDMF( fields, 2, 0, 0, 0 );
	WriteXDMFFooter( 0 );

	delete mesh;
	delete none;
	delete phiClean;
	delete phiModal;
	delete phiNodal;
	delete phiModalVec;
	delete phiNodalVec;
	delete he;
	VecDestroy( nVec );
	VecDestroy( mVec );
	MatDestroy( M2N );

	return EXIT_SUCCESS;
}
