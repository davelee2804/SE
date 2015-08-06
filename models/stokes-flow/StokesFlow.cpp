#include <string>
#include <cmath>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "QuadPoint.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Operator.h"
#include "RHSOp.h"
#include "Vector.h"
#include "Matrix.h"
#include "StokesEqn.h"

using namespace std;
using std::string;

typedef double (AnalyticFunc) (double* coord);

/*
double f_x( double* coord ) { return 2.0 + M_PI*cos(M_PI*coord[0])*sin(M_PI*coord[1]); }
double f_y( double* coord ) { return 0.0 + M_PI*sin(M_PI*coord[0])*cos(M_PI*coord[1]); }
double v_x( double* coord ) { return 1.0 - coord[1]*coord[1]; }
double v_y( double* coord ) { return 0.0; }
double p( double* coord ) { return sin(M_PI*coord[0])*sin(M_PI*coord[1]); }
*/
double v_x( double* coord ) {
	double x = coord[0], y = coord[1];
	double E = (4.0*M_PI*M_PI + 2.0)*exp( 2.0*M_PI ) - exp( 4.0*M_PI ) - 1.0;
	double A = (exp(2.0*M_PI) - 1.0)*exp(M_PI)/E;
	double B = -A;
	double C = (2.0*M_PI - exp(2.0*M_PI) + 1.0)*exp(M_PI)/E;
	double D = -(2.0*M_PI*exp(2.0*M_PI) - exp(2.0*M_PI) + 1.0)*exp(M_PI)/E;

	return sin(M_PI*x)*( ( A*M_PI + C + C*M_PI*y )*exp(M_PI*y) - ( B*M_PI - D + D*M_PI*y )*exp(-M_PI*y) );
}

double v_y( double* coord ) {
	double x = coord[0], y = coord[1];
	double E = (4.0*M_PI*M_PI + 2.0)*exp( 2.0*M_PI ) - exp( 4.0*M_PI ) - 1.0;
	double A = (exp(2.0*M_PI) - 1.0)*exp(M_PI)/E;
	double B = -A;
	double C = (2.0*M_PI - exp(2.0*M_PI) + 1.0)*exp(M_PI)/E;
	double D = -(2.0*M_PI*exp(2.0*M_PI) - exp(2.0*M_PI) + 1.0)*exp(M_PI)/E;

        return -M_PI*cos(M_PI*x)*( (A + C*y)*exp(M_PI*y) + (B + D*y)*exp(-M_PI*y) );
}

double p( double* coord ) {
	double x = coord[0], y = coord[1];
	double E = (4.0*M_PI*M_PI + 2.0)*exp( 2.0*M_PI ) - exp( 4.0*M_PI ) - 1.0;
	double C = (2.0*M_PI - exp(2.0*M_PI) + 1.0)*exp(M_PI)/E;
	double D = -(2.0*M_PI*exp(2.0*M_PI) - exp(2.0*M_PI) + 1.0)*exp(M_PI)/E;

	return -2.0*M_PI*cos(M_PI*x)*( C*exp(M_PI*y) + D*exp(-M_PI*y) );
}

double v_x_top( double* coord ) {
	double min = 0.0, max = 1.0;
	double length = max - min;
	double x = (coord[0] - min)/length;

	return sin(M_PI*x);
}

double IntegrateError( Field* field, AnalyticFunc* func, int dof ) {
	int el_i, pt_i;
	double numeric[2], analytic, gCoord[2], *lCoord, detJac, weight, errorSq, analyticSq;

	errorSq = analyticSq = 0.0;

	for( el_i = 0; el_i < field->mesh->nElsTotal; el_i++ ) {
		for( pt_i = 0; pt_i < field->mesh->el->nPoints; pt_i++ ) {
			lCoord = field->mesh->el->quadPts[pt_i]->coord;
			weight = field->mesh->el->quadPts[pt_i]->weight;
			detJac = field->mesh->DetJacAtCoord( el_i, lCoord );
			field->InterpLocal( el_i, lCoord, numeric );
			field->mesh->LocalToGlobal( lCoord, el_i, gCoord );
			analytic = func( gCoord );
			errorSq += detJac*weight*( analytic - numeric[dof] )*( analytic - numeric[dof] );
			analyticSq += detJac*weight*analytic*analytic;
		}
	}

	return sqrt( errorSq/analyticSq );
}

double IntegrateDivergence( Field* velocity ) {
	int el_i, pt_i;
	double divVel = 0.0, gCoord[2], *lCoord, detJac, weight;
	double** gradV = new double*[2];
	gradV[0] = new double[2];
	gradV[1] = new double[2];

	for( el_i = 0; el_i < velocity->mesh->nElsTotal; el_i++ ) {
		for( pt_i = 0; pt_i < velocity->mesh->el->nPoints; pt_i++ ) {
			lCoord = velocity->mesh->el->quadPts[pt_i]->coord;
			weight = velocity->mesh->el->quadPts[pt_i]->weight;
			detJac = velocity->mesh->DetJacAtCoord( el_i, lCoord );
			velocity->mesh->LocalToGlobal( lCoord, el_i, gCoord );
			velocity->InterpDerivsGlobal( gCoord, gradV );
			divVel += detJac*weight*( gradV[0][0] + gradV[1][1] );
		}
	}

	delete[] gradV[0];
	delete[] gradV[1];
	delete[] gradV;

	return divVel;
}

int main( int argc, char** argv ) {
	int 		N 		= 8;
	int 		nx[2] 		= { 8, 8 };
	double 		min[2] 		= { 0.0, 0.0 }; //{ -1.0, -1.0 };
	double 		max[2] 		= { 1.0, 1.0 };
	bool 		periodic[2] 	= { false, false };
	bool 		vertBCs[2] 	= { true, true };
	bool 		horizBCs[2] 	= { true, false }; //{ true, true };
	Mesh* 		vMesh 		= new Mesh( "vMesh", nx, "legendre", N, min, max, periodic );
	Mesh* 		pMesh 		= new Mesh( "hMesh", nx, "legendre", N-2, min, max, periodic );
	BCs* 		velocityBCs 	= new BCs( vertBCs, vertBCs, horizBCs, horizBCs, vMesh, 0 );
	int 		vSize		= 2*vMesh->nVertsTotal - velocityBCs->size[0] - velocityBCs->size[1];
	BCs* 		pressureBCs	= new BCs( false, false, false, false, pMesh, vSize );
	Field* 		velocity	= new Field( "velocity", vMesh, 2, velocityBCs );
	Field* 		pressure	= new Field( "pressure", pMesh, 1, pressureBCs );
	Field* 		forcing		= new Field( "forcing", vMesh, 2, NULL );
	Field* 		vAnalytic	= new Field( "vAnalytic", vMesh, 2, NULL );
	Field*	 	pAnalytic	= new Field( "pAnalytic", pMesh, 1, NULL );
	StokesEqn* 	stokesEqn;
	Field* 		fields[3];
	char 		tag[] 		= "petsc";

	PetscInitialize( &argc, &argv, (char*)0, tag );

	/* setup */
/*
	velocity->SetBCFunc( "bottom", 0, v_x );
	velocity->SetBCFunc( "top", 0, v_x );
	velocity->SetBCFunc( "left", 0, v_x );
	velocity->SetBCFunc( "right", 0, v_x );
	velocity->SetBCFunc( "bottom", 1, v_y );
	velocity->SetBCFunc( "top", 1, v_y );
	velocity->SetBCFunc( "left", 1, v_y );
	velocity->SetBCFunc( "right", 1, v_y );
*/
	velocity->SetBCConst( "bottom", 0, 0.0 );
	velocity->SetBCConst( "bottom", 1, 0.0 );
	velocity->SetBCConst( "left", 0, 0.0 );
	velocity->SetBCConst( "right", 0, 0.0 );
	velocity->SetBCFunc( "top", 0, v_x_top );
	velocity->SetBCConst( "top", 1, 0.0 );
	vAnalytic->SetICFunc( 0, v_x );
	vAnalytic->SetICFunc( 1, v_y );
	pAnalytic->SetICFunc( 0, p );

	/* solve */
	stokesEqn = new StokesEqn( velocity, pressure, forcing );
	stokesEqn->Solve();

	/* write results */
	vMesh->Save();
	pMesh->Save();
	WriteXDMFHeader( 1 );
	fields[0] = velocity;
	fields[1] = vAnalytic;
	fields[2] = forcing;
	WriteXDMF( fields, 3, 1, 0.0, 0.0 );
	fields[0] = pressure;
	fields[1] = pAnalytic;
	WriteXDMF( fields, 2, 1, 0.0, 0.0 );
	WriteXDMFFooter( 1 );

	cout << "v_x error: " << IntegrateError( velocity, v_x, 0 ) << endl;
	cout << "v_y error: " << IntegrateError( velocity, v_y, 1 ) << endl;
	cout << "p error: " << IntegrateError( pressure, p, 0 ) << endl;
	cout << "net divergence: " << IntegrateDivergence( velocity ) << endl;

	/* cleanup */
	delete stokesEqn;
	delete pAnalytic;
	delete vAnalytic;
	delete forcing;
	delete pressure;
	delete velocity;
	delete pressureBCs;
	delete velocityBCs;
	delete pMesh;
	delete vMesh;

	PetscFinalize();

	return 0;
}
