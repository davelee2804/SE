#include <iostream>
#include <string>
#include <cmath>

#include <fftw3.h>

#include "QuadPoint.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"
#include "Advector.h"
#include "QuasiGeostrophicSpectral.h"

using namespace std;
using std::string;

QuasiGeostrophicSpectral::QuasiGeostrophicSpectral( Field* _vorticity, Field* _streamFunc, Field* _velocity, double _beta ) {
	vorticity  = _vorticity;
	streamFunc = _streamFunc;
	velocity   = _velocity;

	beta = _beta;

	nx = velocity->mesh->nVerts[0] - 1;
	ny = velocity->mesh->nVerts[1] - 1;
	nTot = nx*ny;

	p_real = (double*)fftw_malloc( nTot*sizeof(double) );
	w_real = (double*)fftw_malloc( nTot*sizeof(double) );
	u_real = (double*)fftw_malloc( nTot*sizeof(double) );
	v_real = (double*)fftw_malloc( nTot*sizeof(double) );
	wHat_real = (double*)fftw_malloc( nTot*sizeof(double) );
	p_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	w_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	u_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	v_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	wHat_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );

	p_backward = fftw_plan_dft_c2r_2d( nx, ny, p_fourier, p_real, FFTW_MEASURE );
	w_backward = fftw_plan_dft_c2r_2d( nx, ny, w_fourier, w_real, FFTW_MEASURE );
	u_backward = fftw_plan_dft_c2r_2d( nx, ny, u_fourier, u_real, FFTW_MEASURE );
	v_backward = fftw_plan_dft_c2r_2d( nx, ny, v_fourier, v_real, FFTW_MEASURE );
	wHat_forward = fftw_plan_dft_r2c_2d( nx, ny, wHat_real, wHat_fourier, FFTW_MEASURE );
}

QuasiGeostrophicSpectral::~QuasiGeostrophicSpectral() {
	fftw_free( p_real );
	fftw_free( w_real );
	fftw_free( u_real );
	fftw_free( v_real );
	fftw_free( wHat_real );
	fftw_free( p_fourier );
	fftw_free( w_fourier );
	fftw_free( u_fourier );
	fftw_free( v_fourier );
	fftw_free( wHat_fourier );

	fftw_destroy_plan( p_backward );
	fftw_destroy_plan( w_backward );
	fftw_destroy_plan( u_backward );
	fftw_destroy_plan( v_backward );
	fftw_destroy_plan( wHat_forward );
}

void QuasiGeostrophicSpectral::Solve( double dt, Field* vortPrev, Field* velPrev ) {
	int 		mode_i, node_i;
	int		kx, ky, ksq;
	Advector*	adv;
	double		fourierScale, inv, advFac, nu;

	if( vortPrev && velPrev ) {
		adv = new Advector( vorticity, velocity, vortPrev, velPrev );
		advFac = 1.5;
	}
	else {
		adv = new Advector( vorticity, velocity );
		advFac = 1.0;
	}

	/* 1st order semi-Lagrangian advection in real space */
	adv->cubicInterp = true;
	adv->Advect( dt );
	MapToArray( adv->fieldSL, wHat_real, 0 );
	fftw_execute( wHat_forward );
	fourierScale = 2.0/(nx*ny);
	for( mode_i = 0; mode_i < nx*(ny/2+1); mode_i++ ) {
		wHat_fourier[mode_i][0] *= fourierScale;
		wHat_fourier[mode_i][1] *= fourierScale;
	}

	/* solve the linear terms in fourier space */
	p_fourier[0][0] = p_fourier[0][1] = 0.0;
	w_fourier[0][0] = w_fourier[0][1] = 0.0;
	u_fourier[0][0] = u_fourier[0][1] = 0.0;
	v_fourier[0][0] = v_fourier[0][1] = 0.0;
	for( mode_i = 1; mode_i < nx*(ny/2+1); mode_i++ ) {
                kx = mode_i%(ny/2+1);
                ky = mode_i/(ny/2+1); /* TODO check this! */
                if( ky > ny/2+1 ) {
                        ky = -(ny - ky);
                }
		ksq = kx*kx + ky*ky;
		nu  = 0.01*ksq*ksq;
		inv = 1.0/(kx*kx*dt*dt*beta*beta + advFac*advFac*ksq*ksq - dt*nu);
		p_fourier[mode_i][0] = inv*(  kx*dt*beta*wHat_fourier[mode_i][1] - (advFac*ksq - dt*nu)*wHat_fourier[mode_i][0] );
		p_fourier[mode_i][1] = inv*( -kx*dt*beta*wHat_fourier[mode_i][0] - (advFac*ksq - dt*nu)*wHat_fourier[mode_i][1] );

		w_fourier[mode_i][0] = -ksq*p_fourier[mode_i][0];
		w_fourier[mode_i][1] = -ksq*p_fourier[mode_i][1];

		u_fourier[mode_i][0] =  ky*p_fourier[mode_i][1];
		u_fourier[mode_i][1] = -ky*p_fourier[mode_i][0];

		v_fourier[mode_i][0] = -kx*p_fourier[mode_i][1];
		v_fourier[mode_i][1] =  kx*p_fourier[mode_i][0];
	}

	fftw_execute( p_backward );
	fftw_execute( w_backward );
	fftw_execute( u_backward );
	fftw_execute( v_backward );
	for( node_i = 0; node_i < nTot; node_i++ ) {
		p_real[node_i] *= 0.5;
		w_real[node_i] *= 0.5;
		u_real[node_i] *= 0.5;
		v_real[node_i] *= 0.5;
	}

	MapFromArray( streamFunc, p_real, 0 );
	MapFromArray( vorticity,  w_real, 0 );
	MapFromArray( velocity,   u_real, 0 );
	MapFromArray( velocity,   v_real, 1 );

	delete adv;
}

void QuasiGeostrophicSpectral::MapToArray( Field* field, double* array, int dof ) {
        int node_i, index = 0, topo[2];

        for( node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
                field->mesh->IndexToTopo( node_i, topo );
                if( topo[1] == field->mesh->nVerts[1] - 1 ) {
                        continue;
                }
                if( topo[0] == field->mesh->nVerts[0] - 1 ) {
                        continue;
                }
                array[index] = field->vals[node_i][dof];
                index++;
        }
}

void QuasiGeostrophicSpectral::MapFromArray( Field* field, double* array, int dof ) {
        int node_i, index = 0, topo[2];

        for( node_i = 0; node_i < field->mesh->nVertsTotal; node_i++ ) {
                field->mesh->IndexToTopo( node_i, topo );
                if( topo[1] == field->mesh->nVerts[1] - 1 ) {
                        continue;
                }
                if( topo[0] == field->mesh->nVerts[0] - 1 ) {
                        continue;
                }
                field->vals[node_i][dof] = array[index];
                index++;
        }
        field->PeriodicUpdate();
}
