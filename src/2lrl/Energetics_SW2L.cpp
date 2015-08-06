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
#include "BarotropicEqn.h"
#include "SW2L.h"
#include "Energetics_SW2L.h"

using namespace std;
using std::string;

Energetics_SW2L::Energetics_SW2L( SW2L* _sw, BarotropicEqn* _bt ) {
	sw = _sw;
	bt = _bt;
}

Energetics_SW2L::~Energetics_SW2L() {
}

double Energetics_SW2L::CalcKE( int ts, int layer ) {
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

	bt->CalcHeights( sw->tp, sw->bp, sw->etaInt, &height1, &height2 );

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

double Energetics_SW2L::CalcPE( int ts ) {
	double PE = 0.0, PE_i, e, detJac, weight, *coord;

	sw->etaInt->Read( ts );

	for( int el_i = 0; el_i < sw->etaInt->mesh->nElsTotal; el_i++ ) {
		PE_i = 0.0;
		for( int pt_i = 0; pt_i < sw->etaInt->mesh->el->nPoints; pt_i++ ) {
			weight = sw->etaInt->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->etaInt->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->etaInt->mesh->DetJac( el_i, pt_i );
			sw->etaInt->InterpLocal( el_i, coord, &e );
			PE_i += 0.5*detJac*weight*sw->bp->hFac*e*e;
		}
		PE += PE_i;
	}

	return PE;
}

double Energetics_SW2L::CalcPower_SurfPres( int ts, int layer ) {
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
	bt->pressure->Read( ts );

	bt->CalcHeights( sw->tp, sw->bp, sw->etaInt, &height1, &height2 );

	if( layer == 1 ) {
		velocity = sw->velTop;
		height   = height1;
	}
	else if( layer == 2 ) {
		velocity = sw->velBot;
		height   = height2;
	}

	for( int el_i = 0; el_i < bt->pressure->mesh->nElsTotal; el_i++ ) {
		sp_i = 0.0;
		for( int pt_i = 0; pt_i < bt->pressure->mesh->el->nPoints; pt_i++ ) {
			weight = bt->pressure->mesh->el->quadPts[pt_i]->weight;
			coord  = bt->pressure->mesh->el->quadPts[pt_i]->coord;
			detJac = bt->pressure->mesh->DetJac( el_i, pt_i );
			velocity->InterpLocal( el_i, coord, v );
			height->InterpLocal( el_i, coord, &h );
			bt->pressure->mesh->LocalToGlobal( coord, el_i, gCoord );
			bt->pressure->InterpDerivsGlobal( gCoord, dp );
			sp_i -= detJac*weight*bt->params->pFac*h*( v[0]*dp[0][0] + v[1]*dp[0][1] );
			if( layer == 2 ) {
				sw->etaInt->InterpDerivsGlobal( gCoord, de );
				sp_i -= detJac*weight*sw->bp->hFac*h*( v[0]*de[0][0] + v[1]*de[0][1] );
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

double Energetics_SW2L::CalcPower_WindStress( int ts ) {
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
			tau = sw->tp->tau*cos( sw->tp->k*gCoord[1] );
			ws_i -= detJac*weight*tau*v[0];
		}
		ws += ws_i;
	}

	return ws;
}

double Energetics_SW2L::CalcPower_Friction( int ts ) {
	Field* height1;
	Field* height2;
	double fr = 0.0, fr_i, detJac, weight, *coord, h, v[2];

	sw->velBot->Read( ts );
	sw->etaInt->Read( ts );

	bt->CalcHeights( sw->tp, sw->bp, sw->etaInt, &height1, &height2 );

	for( int el_i = 0; el_i < sw->velBot->mesh->nElsTotal; el_i++ ) {
		fr_i = 0.0;
		for( int pt_i = 0; pt_i < sw->velBot->mesh->el->nPoints; pt_i++ ) {
			weight = sw->velBot->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->velBot->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->velBot->mesh->DetJac( el_i, pt_i );
			sw->velBot->InterpLocal( el_i, coord, v );
			height2->InterpLocal( el_i, coord, &h );
			fr_i -= detJac*weight*sw->bp->gamma*h*( v[0]*v[0] + v[1]*v[1] );
		}
		fr += fr_i;
	}

	delete height1;
	delete height2;

	return fr;
}

double Energetics_SW2L::CalcPower_Viscosity( int ts, int layer ) {
	//TODO: can't take 2nd derivatives of velocity field at present.
	return -1.0;
}
