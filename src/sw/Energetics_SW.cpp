#include <iostream>
#include <string>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>

#include "Base.h"
#include "SWOps.h"
#include "ShallowWaterEqn.h"
#include "Energetics_SW.h"

using namespace std;
using std::string;

Energetics_SW::Energetics_SW( ShallowWaterEqn* _sw ) {
	sw = _sw;
}

Energetics_SW::~Energetics_SW() {}

double Energetics_SW::CalcKE( int ts ) {
	Field* height		= new Field( "sw-height", sw->pressure->mesh, 1, NULL );
	double KE = 0.0, KE_i, v[2], h, detJac, weight, *coord, b = 0.0;

	sw->velocity->Read( ts );
	sw->pressure->Read( ts );

	for( int node_i = 0; node_i < sw->pressure->mesh->nVertsTotal; node_i++ ) {
		if( sw->topography != NULL ) {
			sw->topography->InterpGlobal( sw->pressure->mesh->verts[node_i], &b );
		}
		height->vals[node_i][0] = sw->pressure->vals[node_i][0] + 1.0 - b;
	}

	for( int el_i = 0; el_i < sw->velocity->mesh->nElsTotal; el_i++ ) {
		KE_i = 0.0;
		for( int pt_i = 0; pt_i < sw->velocity->mesh->el->nPoints; pt_i++ ) {
			weight = sw->velocity->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->velocity->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->velocity->mesh->DetJac( el_i, pt_i );
			sw->velocity->InterpLocal( el_i, coord, v );
			height->InterpLocal( el_i, coord, &h );
			KE_i += 0.5*detJac*weight*h*( v[0]*v[0] + v[1]*v[1] );
		}
		KE += KE_i;
	}

	delete height;

	return KE;
}

double Energetics_SW::CalcPE( int ts ) {
	double PE = 0.0, PE_i, e, detJac, weight, *coord;

	sw->pressure->Read( ts );

	for( int el_i = 0; el_i < sw->pressure->mesh->nElsTotal; el_i++ ) {
		PE_i = 0.0;
		for( int pt_i = 0; pt_i < sw->pressure->mesh->el->nPoints; pt_i++ ) {
			weight = sw->pressure->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->pressure->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->pressure->mesh->DetJac( el_i, pt_i );
			sw->pressure->InterpLocal( el_i, coord, &e );
			PE_i += 0.5*detJac*weight*sw->params->g*e*e;
		}
		PE += PE_i;
	}

	return PE;
}

double Energetics_SW::CalcPower_WindStress( int ts ) {
	double	ws = 0.0, ws_i, detJac, weight, *coord, v[2], gCoord[2], tau;

	sw->velocity->Read( ts );

	for( int el_i = 0; el_i < sw->velocity->mesh->nElsTotal; el_i++ ) {
		ws_i = 0.0;
		for( int pt_i = 0; pt_i < sw->velocity->mesh->el->nPoints; pt_i++ ) {
			weight = sw->velocity->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->velocity->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->velocity->mesh->DetJac( el_i, pt_i );
			sw->velocity->InterpLocal( el_i, coord, v );
			sw->velocity->mesh->LocalToGlobal( coord, el_i, gCoord );
			tau = sw->params->tau*cos( sw->params->kws*gCoord[1] );
			ws_i -= detJac*weight*tau*v[0];
		}
		ws += ws_i;
	}

	return ws;
}

double Energetics_SW::CalcPower_Friction( int ts ) {
	Field* height		= new Field( "sw-height", sw->pressure->mesh, 1, NULL );
	double fr = 0.0, fr_i, detJac, weight, *coord, h, v[2], b = 0.0;

	sw->velocity->Read( ts );
	sw->pressure->Read( ts );

	for( int node_i = 0; node_i < sw->pressure->mesh->nVertsTotal; node_i++ ) {
		if( sw->topography != NULL ) {
			sw->topography->InterpGlobal( sw->pressure->mesh->verts[node_i], &b );
		}
		height->vals[node_i][0] = sw->pressure->vals[node_i][0] + 1.0 - b;
	}

	for( int el_i = 0; el_i < sw->velocity->mesh->nElsTotal; el_i++ ) {
		fr_i = 0.0;
		for( int pt_i = 0; pt_i < sw->velocity->mesh->el->nPoints; pt_i++ ) {
			weight = sw->velocity->mesh->el->quadPts[pt_i]->weight;
			coord  = sw->velocity->mesh->el->quadPts[pt_i]->coord;
			detJac = sw->velocity->mesh->DetJac( el_i, pt_i );
			sw->velocity->InterpLocal( el_i, coord, v );
			height->InterpLocal( el_i, coord, &h );
			fr_i -= detJac*weight*sw->params->gamma*h*( v[0]*v[0] + v[1]*v[1] );
		}
		fr += fr_i;
	}

	delete height;

	return fr;
}

double Energetics_SW::CalcPower_Viscosity( int ts ) {
	//TODO: can't take 2nd derivatives of velocity field at present.
	return -1.0;
}
