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
#include "NavierStokesSpectral.h"

using namespace std;
using std::string;

NavierStokesSpectral::NavierStokesSpectral( Field* _velocity, Field* _pressure, double _nu ) {
	velocity = _velocity;
	pressure = _pressure;
	nu = _nu;

	nx = velocity->mesh->nVerts[0] - 1;
	ny = velocity->mesh->nVerts[1] - 1;
	nTot = nx*ny;

	u_real = (double*)fftw_malloc( nTot*sizeof(double) );
	v_real = (double*)fftw_malloc( nTot*sizeof(double) );
	p_real = (double*)fftw_malloc( nTot*sizeof(double) );
	uHat_real = (double*)fftw_malloc( nTot*sizeof(double) );
	vHat_real = (double*)fftw_malloc( nTot*sizeof(double) );
	u_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	v_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	p_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	uHat_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	vHat_fourier = (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );

	u_backward = fftw_plan_dft_c2r_2d( nx, ny, u_fourier, u_real, FFTW_MEASURE );
	v_backward = fftw_plan_dft_c2r_2d( nx, ny, v_fourier, v_real, FFTW_MEASURE );
	p_backward = fftw_plan_dft_c2r_2d( nx, ny, p_fourier, p_real, FFTW_MEASURE );
	uHat_forward  = fftw_plan_dft_r2c_2d( nx, ny, uHat_real, uHat_fourier, FFTW_MEASURE );
	vHat_forward  = fftw_plan_dft_r2c_2d( nx, ny, vHat_real, vHat_fourier, FFTW_MEASURE );
}

NavierStokesSpectral::~NavierStokesSpectral() {
	fftw_free( u_real );
	fftw_free( v_real );
	fftw_free( p_real );
	fftw_free( uHat_real );
	fftw_free( vHat_real );
	fftw_free( u_fourier );
	fftw_free( v_fourier );
	fftw_free( p_fourier );
	fftw_free( uHat_fourier );
	fftw_free( vHat_fourier );

	fftw_destroy_plan( u_backward );
	fftw_destroy_plan( v_backward );
	fftw_destroy_plan( p_backward );
	fftw_destroy_plan( uHat_forward );
	fftw_destroy_plan( vHat_forward );
}

void NavierStokesSpectral::Solve( double dt, Field* velPrev ) {
	if( velPrev ) {
		SolveSecondOrder( dt, velPrev );
	}
	else {
		SolveFirstOrder( dt );
	}
}

void NavierStokesSpectral::SolveFirstOrder( double dt ) {
	int 		mode_i, node_i;
	int		kx, ky, ksq;
	Advector*	adv	= new Advector( velocity, velocity );
	fftw_complex*	uHatHat	= (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	fftw_complex*	vHatHat	= (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	double		eta, a, b, c, detInv;
	double		fourierScale;

	/* 1st order semi-Lagrangian advection in real space */
	adv->Advect( dt );
	MapToArray( adv->fieldSL, uHat_real, 0 );
	MapToArray( adv->fieldSL, vHat_real, 1 );
	fftw_execute( uHat_forward );
	fftw_execute( vHat_forward );
	fourierScale = 2.0/(nx*ny);
	for( mode_i = 0; mode_i < nx*(ny/2+1); mode_i++ ) {
		uHat_fourier[mode_i][0] *= fourierScale;
		uHat_fourier[mode_i][1] *= fourierScale;
		vHat_fourier[mode_i][0] *= fourierScale;
		vHat_fourier[mode_i][1] *= fourierScale;
	}

	/* solve the pressure poisson eqn in fourier space */
	p_fourier[0][0] = p_fourier[0][1] = 0.0;
	for( mode_i = 1; mode_i < nx*(ny/2+1); mode_i++ ) {
                kx = mode_i%(ny/2+1);
                ky = mode_i/(ny/2+1); /* TODO check this! */
                if( ky > ny/2+1 ) {
                        ky = -(ny - ky);
                }
		ksq = kx*kx + ky*ky;
		p_fourier[mode_i][0] = +(kx*uHat_fourier[mode_i][1] + ky*vHat_fourier[mode_i][1])/(dt*ksq);
		p_fourier[mode_i][1] = -(kx*uHat_fourier[mode_i][0] + ky*vHat_fourier[mode_i][0])/(dt*ksq);
	}

	/* velocity update */
	for( mode_i = 1; mode_i < nx*(ny/2+1); mode_i++ ) {
                kx = mode_i%(ny/2+1);
                ky = mode_i/(ny/2+1); /* TODO check this! */
                if( ky > ny/2+1 ) {
                        ky = -(ny - ky);
                }
		uHatHat[mode_i][0] = uHat_fourier[mode_i][0] + dt*kx*p_fourier[mode_i][1];
		uHatHat[mode_i][1] = uHat_fourier[mode_i][1] - dt*kx*p_fourier[mode_i][0];
		vHatHat[mode_i][0] = vHat_fourier[mode_i][0] + dt*ky*p_fourier[mode_i][1];
		vHatHat[mode_i][1] = vHat_fourier[mode_i][1] - dt*ky*p_fourier[mode_i][0];
	}

	/* solve the velocity helmholtz eqn in fourier space */
	u_fourier[0][0] = u_fourier[0][1] = v_fourier[0][0] = v_fourier[0][1] = 0.0;
	for( mode_i = 1; mode_i < nx*(ny/2+1); mode_i++ ) {
                kx = mode_i%(ny/2+1);
                ky = mode_i/(ny/2+1); /* TODO check this! */
                if( ky > ny/2+1 ) {
                        ky = -(ny - ky);
                }
		ksq = kx*kx + ky*ky;

		eta = nu + SVV( kx, ky, (7*nx)/8 );
		a = 1.0 + dt*eta*(2*kx*kx + ky*ky);
		b = dt*eta*kx*ky;
		c = 1.0 + dt*eta*(2*ky*ky + kx*kx);
		detInv = 1.0/(a*c - b*b);
		u_fourier[mode_i][0] = detInv*(c*uHatHat[mode_i][0] - b*vHatHat[mode_i][0]);
		u_fourier[mode_i][1] = detInv*(c*uHatHat[mode_i][1] - b*vHatHat[mode_i][1]);
		v_fourier[mode_i][0] = detInv*(a*vHatHat[mode_i][0] - b*uHatHat[mode_i][0]);
		v_fourier[mode_i][1] = detInv*(a*vHatHat[mode_i][1] - b*uHatHat[mode_i][1]);
	}

	fftw_execute( u_backward );
	fftw_execute( v_backward );
	fftw_execute( p_backward );
	for( node_i = 0; node_i < nTot; node_i++ ) {
		u_real[node_i] *= 0.5;
		v_real[node_i] *= 0.5;
		p_real[node_i] *= 0.5;
	}

	MapFromArray( velocity, u_real, 0 );
	MapFromArray( velocity, v_real, 1 );
	MapFromArray( pressure, p_real, 0 );

	delete adv;
	fftw_free( uHatHat );
	fftw_free( vHatHat );
}

void NavierStokesSpectral::SolveSecondOrder( double dt, Field* velPrev ) {
	int 		mode_i, node_i;
	int		kx, ky, ksq;
	Advector*	adv	= new Advector( velocity, velocity, velPrev, velPrev );
	Field*		velHat	= new Field( "velHat", velocity->mesh, 2, NULL );
	fftw_complex*	uHatHat	= (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	fftw_complex*	vHatHat	= (fftw_complex*)fftw_malloc( nx*(ny/2+1)*sizeof(fftw_complex) );
	double		eta, a, b, c, detInv;
	double		fourierScale;

	/* 1st order semi-Lagrangian advection in real space */
	adv->Advect( dt );
	for( node_i = 0; node_i < velocity->mesh->nVertsTotal; node_i++ ) {
		velHat->vals[node_i][0] = adv->fieldSL->vals[node_i][0];
		velHat->vals[node_i][1] = adv->fieldSL->vals[node_i][1];
	}
	MapToArray( velHat, uHat_real, 0 );
	MapToArray( velHat, vHat_real, 1 );
	fftw_execute( uHat_forward );
	fftw_execute( vHat_forward );
	fourierScale = 2.0/(nx*ny);
	for( mode_i = 0; mode_i < nx*(ny/2+1); mode_i++ ) {
		uHat_fourier[mode_i][0] *= fourierScale;
		uHat_fourier[mode_i][1] *= fourierScale;
		vHat_fourier[mode_i][0] *= fourierScale;
		vHat_fourier[mode_i][1] *= fourierScale;
	}

	/* solve the pressure poisson eqn in fourier space */
	p_fourier[0][0] = p_fourier[0][1] = 0.0;
	for( mode_i = 1; mode_i < nx*(ny/2+1); mode_i++ ) {
                kx = mode_i%(ny/2+1);
                ky = mode_i/(ny/2+1); /* TODO check this! */
                if( ky > ny/2+1 ) {
                        ky = -(ny - ky);
                }
		ksq = kx*kx + ky*ky;
		p_fourier[mode_i][0] = +(kx*uHat_fourier[mode_i][1] + ky*vHat_fourier[mode_i][1])/(dt*ksq);
		p_fourier[mode_i][1] = -(kx*uHat_fourier[mode_i][0] + ky*vHat_fourier[mode_i][0])/(dt*ksq);
	}

	/* velocity update */
	for( mode_i = 1; mode_i < nx*(ny/2+1); mode_i++ ) {
                kx = mode_i%(ny/2+1);
                ky = mode_i/(ny/2+1); /* TODO check this! */
                if( ky > ny/2+1 ) {
                        ky = -(ny - ky);
                }
		uHatHat[mode_i][0] = uHat_fourier[mode_i][0] + dt*kx*p_fourier[mode_i][1];
		uHatHat[mode_i][1] = uHat_fourier[mode_i][1] - dt*kx*p_fourier[mode_i][0];
		vHatHat[mode_i][0] = vHat_fourier[mode_i][0] + dt*ky*p_fourier[mode_i][1];
		vHatHat[mode_i][1] = vHat_fourier[mode_i][1] - dt*ky*p_fourier[mode_i][0];
	}

	/* solve the velocity helmholtz eqn in fourier space */
	u_fourier[0][0] = u_fourier[0][1] = v_fourier[0][0] = v_fourier[0][1] = 0.0;
	for( mode_i = 1; mode_i < nx*(ny/2+1); mode_i++ ) {
                kx = mode_i%(ny/2+1);
                ky = mode_i/(ny/2+1); /* TODO check this! */
                if( ky > ny/2+1 ) {
                        ky = -(ny - ky);
                }
		ksq = kx*kx + ky*ky;

		eta = nu + SVV( kx, ky, (7*nx)/8 );
		a = 1.5 + dt*eta*(2*kx*kx + ky*ky);
		b = dt*eta*kx*ky;
		c = 1.5 + dt*eta*(2*ky*ky + kx*kx);
		detInv = 1.0/(a*c - b*b);
		u_fourier[mode_i][0] = detInv*(c*uHatHat[mode_i][0] - b*vHatHat[mode_i][0]);
		u_fourier[mode_i][1] = detInv*(c*uHatHat[mode_i][1] - b*vHatHat[mode_i][1]);
		v_fourier[mode_i][0] = detInv*(a*vHatHat[mode_i][0] - b*uHatHat[mode_i][0]);
		v_fourier[mode_i][1] = detInv*(a*vHatHat[mode_i][1] - b*uHatHat[mode_i][1]);
	}

	fftw_execute( u_backward );
	fftw_execute( v_backward );
	fftw_execute( p_backward );
	for( node_i = 0; node_i < nTot; node_i++ ) {
		u_real[node_i] *= 0.5;
		v_real[node_i] *= 0.5;
		p_real[node_i] *= 0.5;
	}

	MapFromArray( velocity, u_real, 0 );
	MapFromArray( velocity, v_real, 1 );
	MapFromArray( pressure, p_real, 0 );

	delete adv;
	delete velHat;
	fftw_free( uHatHat );
	fftw_free( vHatHat );
}

void NavierStokesSpectral::MapToArray( Field* field, double* array, int dof ) {
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

void NavierStokesSpectral::MapFromArray( Field* field, double* array, int dof ) {
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

/* assumes ns = ny! */
double NavierStokesSpectral::SVV( int kx, int ky, int cutoff ) {
	double	eps, alpha;

/*
	if( kx <= cutoff && ky <= cutoff ) {
		return 0.0;
	}

	eps = 1.0/nx;
	alpha = ((double)((kx - nx)*(ky - ny)))/((kx - cutoff)*(ky - cutoff));
	return eps*exp(-alpha);
*/
	if( kx + ky <= cutoff ) {
		return 0.0;
	}

	eps = 1.0/nx;
	alpha = ((double)(kx + ky - nx))/(kx + ky - cutoff);
	return eps*exp(-alpha*alpha);
}
