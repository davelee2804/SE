#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "LinSys.h"
#include "NavierStokesEqn.h"

using namespace std;
using std::string;

/* solves the 2D incompressible Navier-Stokes equations
   reference:
	Xu, J. D. Xiu & G. E. Karniadakis (2002) "A Semi-Lagrangian Method for Turbulence 
	Simulations Using Mixed Spectral Discretisations" Journal of Scientific Computing
	vol. 17, 585 - 597 */

NavierStokesEqn::NavierStokesEqn( Field* velocity, Field* pressure, double _dt ) {
	dt = _dt;
	secondOrderDone = false;

	velStep2 = new Field( "velStep2", velocity->mesh, 2, NULL );
	divVelStep1 = new Field( "divVelStep1", pressure->mesh, 1, NULL );

	poisson   = new PoissonEqn( pressure, divVelStep1, -dt );
	helmholtz = NULL;
}

NavierStokesEqn::~NavierStokesEqn() {
	delete velStep2;
	delete divVelStep1;
	delete poisson;
	delete helmholtz;
}

void NavierStokesEqn::Solve( Field* velocity, Field* pressure, Field* prevVel, double nu, int svvCutoff ) {
	if( prevVel ) {
		SolveSecondOrder( velocity, pressure, prevVel, nu, svvCutoff );
	}
	else {
		SolveFirstOrder( velocity, pressure, nu, svvCutoff );
	}
}

void NavierStokesEqn::SolveFirstOrder( Field* velocity, Field* pressure, double nu, int svvCutoff ) {
	Advector* 	adv;
	int		node_i;
	double		**dp, **dv;

	firstOrderInit = true;

	cout << "advective step...\n";
	adv = new Advector( velocity, velocity );
	adv->Advect( dt );

	cout << "pressure step...\n";

	dv = new double*[2];
	dv[0] = new double[2];
	dv[1] = new double[2];
	for( node_i = 0; node_i < pressure->mesh->nVertsTotal; node_i++ ) {
		adv->fieldSL->InterpDerivsGlobal( pressure->mesh->verts[node_i], dv );
		divVelStep1->vals[node_i][0] = dv[0][0] + dv[1][1];
	}
	delete[] dv[0];
	delete[] dv[1];
	delete[] dv;

	/* remove the null space */
	poisson->Solve( true );

	dp = new double*[1];
	dp[0] = new double[2];
	for( node_i = 0; node_i < velStep2->mesh->nVertsTotal; node_i++ ) {
		pressure->InterpDerivsGlobal( velStep2->mesh->verts[node_i], dp );
		velStep2->vals[node_i][0] = adv->fieldSL->vals[node_i][0] - dt*dp[0][0];
		velStep2->vals[node_i][1] = adv->fieldSL->vals[node_i][1] - dt*dp[0][1];
	}
	delete[] dp[0];
	delete[] dp;

	cout << "viscous step...\n";

	if( !helmholtz ) {
		helmholtz = new HelmholtzEqn( velocity, velStep2, 1.0, nu, dt, svvCutoff );
	}
	//helmholtz->Solve( "helm_" );
	helmholtz->Solve( "helm_", false );

	delete adv;
}

void NavierStokesEqn::SolveSecondOrder( Field* velocity, Field* pressure, Field* prevVel, double nu, int svvCutoff ) {
	Advector* 	adv;
	int		node_i;
	double		**dp, **dv;

	cout << "advective step...\n";
	adv = new Advector( velocity, velocity, prevVel, prevVel );
	adv->Advect( dt );

	cout << "pressure step...\n";

	dv = new double*[2];
	dv[0] = new double[2];
	dv[1] = new double[2];
	for( node_i = 0; node_i < pressure->mesh->nVertsTotal; node_i++ ) {
		adv->fieldSL->InterpDerivsGlobal( pressure->mesh->verts[node_i], dv );
		divVelStep1->vals[node_i][0] = (dv[0][0] + dv[1][1])*1.5*10.0;
	}
	delete[] dv[0];
	delete[] dv[1];
	delete[] dv;

	/* remove the null space */
	poisson->Solve( true );

	dp = new double*[1];
	dp[0] = new double[2];
	for( node_i = 0; node_i < velStep2->mesh->nVertsTotal; node_i++ ) {
		pressure->InterpDerivsGlobal( velStep2->mesh->verts[node_i], dp );
		velStep2->vals[node_i][0] = adv->fieldSL->vals[node_i][0] - dt*dp[0][0];
		velStep2->vals[node_i][1] = adv->fieldSL->vals[node_i][1] - dt*dp[0][1];
	}
	delete[] dp[0];
	delete[] dp;

	cout << "viscous step...\n";

	if( !secondOrderDone ) {
		if( helmholtz ) {
			delete helmholtz;
		}
		helmholtz = new HelmholtzEqn( velocity, velStep2, 1.5, nu, dt, svvCutoff );
	}
	//helmholtz->Solve( "helm_" );
	helmholtz->Solve( "helm_", false );

	secondOrderDone = true;

	delete adv;
}
