#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SW.h"
#include "Params.h"
#include "2LRLOps.h"
#include "TwoLayerSW.h"
#include "TwoLayerSW_Energetics.h"

using namespace std;
using std::string;

TwoLayerSW_Energetics::TwoLayerSW_Energetics( TwoLayerSW* _sw ) {
	sw = _sw;
}

TwoLayerSW_Energetics::~TwoLayerSW_Energetics() {}

double TwoLayerSW_Energetics::CalcKE( int ts, int layer ) {
	Field*	velocity	= NULL;
	Field*	height		= NULL;
	Field*	height1;
	Field*	height2;
	double	KE = 0.0, KE_i, v[2], h, detJac, weight, *coord;

	if( layer != 1 && layer != 2 ) {
		return -1.0;
	}

	sw->velTop->Read( ts );
	sw->velBot->Read( ts );
	sw->etaInt->Read( ts );

	CalcHeights( sw->params, sw->etaInt, &height1, &height2 );

	if( layer == 1 ) {
		velocity = sw->velTop;
		height   = height1;
	}
	else if( layer == 2 ) {
		velocity = sw->velBot;
		height   = height2;
	}

	for( int el_i = 0; el_i < velocity->mesh->nElsTotal; el_i++ ) {
		KE_i = 0.0;
		for( int pt_i = 0; pt_i < velocity->mesh->el->nPoints; pt_i++ ) {
			weight = velocity->mesh->el->quadPts[pt_i]->weight;
			coord  = velocity->mesh->el->quadPts[pt_i]->coord;
			detJac = velocity->mesh->DetJac( el_i, pt_i );
			velocity->InterpLocal( el_i, coord, v );
			height->InterpLocal( el_i, coord, &h );
			KE_i += 0.5*detJac*weight*h*( v[0]*v[0] + v[1]*v[1] );
		}
		KE += KE_i;
	}

	delete height1;
	delete height2;

	return KE;
}

double TwoLayerSW_Energetics::CalcPE( int ts ) {
	double PE = 0.0, PE_i, e, detJac, weight, *coord;

	sw->etaInt->Read( ts );

	for( int el_i = 0; el_i < sw->etaInt->mesh->nElsTotal; el_i++ ) {
		PE_i = 0.0;
		for( int pt_i = 0; pt_i < sw->etaInt->mesh->el->nPoints; pt_i++ ) {
			weight = sw->etaInt->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->etaInt->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->etaInt->mesh->DetJac( el_i, pt_i );
			sw->etaInt->InterpLocal( el_i, coord, &e );
			PE_i += 0.5*detJac*weight*sw->params->g*e*e;
		}
		PE += PE_i;
	}

	return PE;
}

double TwoLayerSW_Energetics::CalcPower_SurfPres( int ts, int layer ) {
	Field*	velocity	= NULL;
	Field*	height		= NULL;
	Field*	height1;
	Field*	height2;
	double	sp = 0.0, sp_i, **dp, **de, h, v[2], detJac, weight, *coord, gCoord[2];

	if( layer != 1 && layer != 2 ) {
		return -1.0;
	}

	dp = new double*[1];
	dp[0] = new double[2];
	de = new double*[1];
	de[0] = new double[2];

	sw->velTop->Read( ts );
	sw->velBot->Read( ts );
	sw->etaInt->Read( ts );
	sw->presSurf->Read( ts );

	CalcHeights( sw->params, sw->etaInt, &height1, &height2 );

	if( layer == 1 ) {
		velocity = sw->velTop;
		height   = height1;
	}
	else if( layer == 2 ) {
		velocity = sw->velBot;
		height   = height2;
	}

	for( int el_i = 0; el_i < sw->presSurf->mesh->nElsTotal; el_i++ ) {
		sp_i = 0.0;
		for( int pt_i = 0; pt_i < sw->presSurf->mesh->el->nPoints; pt_i++ ) {
			weight = sw->presSurf->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->presSurf->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->presSurf->mesh->DetJac( el_i, pt_i );
			velocity->InterpLocal( el_i, coord, v );
			height->InterpLocal( el_i, coord, &h );
			sw->presSurf->mesh->LocalToGlobal( coord, el_i, gCoord );
			sw->presSurf->InterpDerivsGlobal( gCoord, dp );
			sp_i -= detJac*weight*sw->params->p*h*( v[0]*dp[0][0] + v[1]*dp[0][1] );
			if( layer == 2 ) {
				sw->etaInt->InterpDerivsGlobal( gCoord, de );
				sp_i -= detJac*weight*sw->params->H2*h*( v[0]*de[0][0] + v[1]*de[0][1] );
			}
		}
		sp += sp_i;
	}

	delete[] dp[0];
	delete[] dp;
	delete[] de[0];
	delete[] de;
	delete height1;
	delete height2;

	return sp;
}

double TwoLayerSW_Energetics::CalcPower_WindStress( int ts ) {
	double	ws = 0.0, ws_i, detJac, weight, *coord, v[2], gCoord[2], tau;

	sw->velTop->Read( ts );

	for( int el_i = 0; el_i < sw->velTop->mesh->nElsTotal; el_i++ ) {
		ws_i = 0.0;
		for( int pt_i = 0; pt_i < sw->velTop->mesh->el->nPoints; pt_i++ ) {
			weight = sw->velTop->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->velTop->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->velTop->mesh->DetJac( el_i, pt_i );
			sw->velTop->InterpLocal( el_i, coord, v );
			sw->velTop->mesh->LocalToGlobal( coord, el_i, gCoord );
			tau = sw->params->tau*sin( sw->params->k*gCoord[1] );
			ws_i += detJac*weight*tau*v[0];
		}
		ws += ws_i;
	}

	return ws;
}

double TwoLayerSW_Energetics::CalcPower_Friction( int ts ) {
/*
	Field* height1;
	Field* height2;
	double fr = 0.0, fr_i, detJac, weight, *coord, h, v[2];

	sw->velBot->Read( ts );
	sw->etaInt->Read( ts );

	CalcHeights( sw->tp, sw->bp, sw->etaInt, &height1, &height2 );

	for( int el_i = 0; el_i < sw->velBot->mesh->nElsTotal; el_i++ ) {
		fr_i = 0.0;
		for( int pt_i = 0; pt_i < sw->velBot->mesh->el->nPoints; pt_i++ ) {
			weight = sw->velBot->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->velBot->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->velBot->mesh->DetJac( el_i, pt_i );
			sw->velBot->InterpLocal( el_i, coord, v );
			height2->InterpLocal( el_i, coord, &h );
			fr_i -= detJac*weight*sw->params->gamma*h*( v[0]*v[0] + v[1]*v[1] );
		}
		fr += fr_i;
	}

	delete height1;
	delete height2;

	return fr;
*/
	return 0.0;
}

double TwoLayerSW_Energetics::CalcPower_Viscosity( int ts, int layer ) {
	Field* 		height1;
	Field* 		height2;
	Field* 		height;
	Field*		velocity;
	double		**du, u[2], h, visc = 0.0, visc_e, weight, *coord, detJac;

	if( layer != 1 && layer != 2 ) {
		return -1.0;
	}

	sw->velTop->Read( ts );
	sw->velBot->Read( ts );
	sw->etaInt->Read( ts );

	CalcHeights( sw->params, sw->etaInt, &height1, &height2 );

	if( layer == 1 ) {
		velocity = sw->velTop;
		height   = height1;
	}
	else if( layer == 2 ) {
		velocity = sw->velBot;
		height   = height2;
	}

	du = new double*[2];
	du[0] = new double[3];
	du[1] = new double[3];

	for( int el_i = 0; el_i < velocity->mesh->nElsTotal; el_i++ ) {
		visc_e = 0.0;
		for( int pt_i = 0; pt_i < velocity->mesh->el->nNodes; pt_i++ ) {
			weight = velocity->mesh->el->quadPts[pt_i]->weight;
			coord  = velocity->mesh->el->quadPts[pt_i]->coord;
			detJac = velocity->mesh->DetJac( el_i, pt_i );
			velocity->InterpLocal( el_i, coord, u );
			velocity->InterpSecondDerivs( el_i, pt_i, du );
			height->InterpLocal( el_i, coord, &h );
			visc_e += detJac*weight*sw->params->nu*h*( u[0]*(du[0][0]+du[0][1]) + u[1]*(du[1][0]+du[1][1]) );
		}
		visc += visc_e;
	}

	delete height1;
	delete height2;
	delete[] du[0];
	delete[] du[1];
	delete[] du;

	return visc;
}

void TwoLayerSW_Energetics::CalcHeights( TwoLayerParams* params, Field* eta, Field** height1, Field** height2 ) {
	double e;

	*height1 = new Field( "height1", eta->mesh, 1, NULL );
	*height2 = new Field( "height2", eta->mesh, 1, NULL );

	for( int node_i = 0; node_i < eta->mesh->nVertsTotal; node_i++ ) {
		e = eta->vals[node_i][0];
		(*height1)->vals[node_i][0] = params->H1 - e;
		(*height2)->vals[node_i][0] = params->H2 + e;
	}
}
