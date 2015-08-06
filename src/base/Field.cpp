#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <hdf5.h>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"

using namespace std;
using std::string;

Field::Field( string _name, Mesh* _mesh, int _nDofs, BCs* _bcs ) {
	name = _name;
	mesh = _mesh;
	nDofs = _nDofs;
	bcs = _bcs;

	vals = new double*[mesh->nVertsTotal];
	for( int vert_i = 0; vert_i < mesh->nVertsTotal; vert_i++ ) {
		vals[vert_i] = new double[nDofs];
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			vals[vert_i][dof_i] = 0.0;
		}
	}

	GNix = new double*[mesh->dim];
	for( int dim_i = 0; dim_i < mesh->dim; dim_i++ ) {
		GNix[dim_i] = new double[mesh->el->nNodes];
	}
}

Field::~Field() {
	for( int vert_i = 0; vert_i < mesh->nVertsTotal; vert_i++ ) {
		delete[] vals[vert_i];
	}
	delete[] vals;

	for( int dim_i = 0; dim_i < mesh->dim; dim_i++ ) {
		delete[] GNix[dim_i];
	}
	delete[] GNix;
}

void Field::Copy( Field* field ) {
	for( int node_i = 0; node_i < mesh->nVertsTotal; node_i++ ) {
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			vals[node_i][dof_i] = field->vals[node_i][dof_i];
		}
	}
}

void Field::InterpLocal( int el_i, double* lCoord, double* val ) {
	double* Ni 		= mesh->el->ShapeFuncsAtCoord( lCoord );
	int* 	elNodes 	= mesh->ElNodes( el_i );

	for( int dof_i = 0; dof_i < nDofs; dof_i++ ) { val[dof_i] = 0.0; }

	for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			val[dof_i] += Ni[node_i]*vals[elNodes[node_i]][dof_i];
		}
	}
}

void Field::InterpGlobal( double* gCoord, double* val ) {
	int 	el_i;
	double	lCoord[3];

	mesh->GlobalToLocal( gCoord, &el_i, lCoord );
	InterpLocal( el_i, lCoord, val );
}

void Field::InterpDerivsLocal( int el_i, double* lCoord, double** gradVal ) {
	double**	GNix		= mesh->el->ShapeFuncDerivsAtCoord( lCoord );
	int*		elNodes		= mesh->ElNodes( el_i );

	for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {	
		for( int dim_i = 0; dim_i < mesh->dim; dim_i++ ) {
			gradVal[dof_i][dim_i] = 0.0;
		}
	}
	for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			for( int dim_i = 0; dim_i < mesh->dim; dim_i++ ) {
				gradVal[dof_i][dim_i] += GNix[dim_i][node_i]*vals[elNodes[node_i]][dof_i];
			}
		}
	}
}

void Field::InterpDerivsGlobal( double* gCoord, double** gradVal ) {
	int		el_i;
	double 		lCoord[2], detJac;
	int*		elNodes;
	double**	GNix;
	
	mesh->GlobalToLocal( gCoord, &el_i, lCoord );
	elNodes = mesh->ElNodes( el_i );
	GNix = mesh->ShapeFuncDerivsAtCoord( el_i, lCoord, &detJac );

	for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
		for( int dim_i = 0; dim_i < mesh->dim; dim_i++ ) {
			gradVal[dof_i][dim_i] = 0.0;
		}
	}
	for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			for( int dim_i = 0; dim_i < mesh->dim; dim_i++ ) {
				gradVal[dof_i][dim_i] += GNix[dim_i][node_i]*vals[elNodes[node_i]][dof_i];
			}
		}
	}
}

void Field::InterpSecondDerivs( int el_i, int pt_i, double** g2Val ) {
	int*	elNodes;
	double	detJac, **GNixx;

	elNodes = mesh->ElNodes( el_i );
	GNixx = mesh->ShapeFuncSecondDerivs( el_i, pt_i, &detJac );

	for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
		/* dims are phi,xx, phi,yy, phi,xy */
		for( int dim_i = 0; dim_i < 3; dim_i++ ) {
			g2Val[dof_i][dim_i] = 0.0;
		}
	}
	for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			for( int dim_i = 0; dim_i < 3; dim_i++ ) {
				g2Val[dof_i][dim_i] += GNixx[dim_i][node_i]*vals[elNodes[node_i]][dof_i];
			}
		}
	}
}

void Field::InterpSecondDerivsAtCoord( int el_i, double* coord, double** g2Val ) {
	int*	elNodes;
	double	detJac, **GNixx;

	elNodes = mesh->ElNodes( el_i );
	GNixx = mesh->ShapeFuncSecondDerivsAtCoord( el_i, coord, &detJac );

	for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
		/* dims are phi,xx, phi,yy, phi,xy */
		for( int dim_i = 0; dim_i < 3; dim_i++ ) {
			g2Val[dof_i][dim_i] = 0.0;
		}
	}
	for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			for( int dim_i = 0; dim_i < 3; dim_i++ ) {
				g2Val[dof_i][dim_i] += GNixx[dim_i][node_i]*vals[elNodes[node_i]][dof_i];
			}
		}
	}
}

void Field::InterpDerivsWithGlobalShapeFuncDerivs( int el_i, double** globalGNix, double** gradVal ) {
	int*	elNodes = mesh->ElNodes( el_i );

	for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
		for( int dim_i = 0; dim_i < mesh->dim; dim_i++ ) {
			gradVal[dof_i][dim_i] = 0.0;
		}
	}
	for( int node_i = 0; node_i < mesh->el->nNodes; node_i++ ) {
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			for( int dim_i = 0; dim_i < mesh->dim; dim_i++ ) {
				gradVal[dof_i][dim_i] += globalGNix[dim_i][node_i]*vals[elNodes[node_i]][dof_i];
			}
		}
	}
}

double Field::Integrate( int dof, bool average ) {
        double detJac, weight, x, xtot = 0.0, A1, A2 = 0.0;
        int* elNodes;

        A1 = (mesh->max[0] - mesh->min[0])*(mesh->max[1] - mesh->min[1]);

        for( int el_i = 0; el_i < mesh->nElsTotal; el_i++ ) {
                elNodes = mesh->ElNodes( el_i );
                for( int pt_i = 0; pt_i < mesh->el->nPoints; pt_i++ ) {
                        weight = mesh->el->quadPts[pt_i]->weight;
                        detJac = mesh->DetJac( el_i, pt_i );
                        x      = vals[elNodes[pt_i]][dof];
                        xtot   += detJac*weight*x;
                        A2     += detJac*weight;
                }
        }

        if( average ) { xtot /= A1; }
        return xtot;
}

void Field::SetBCConst( string side, int dof, double val ) {
	int 	node_i, numNodes;
	int* 	nodes 	= bcs->GetSide( side, &numNodes );

	for( node_i = 0; node_i < numNodes; node_i++ ) {
		vals[nodes[node_i]][dof] = val;
	}
}

void Field::SetBCFunc( string side, int dof, BCFunc* bcFunc ) {
	int 	node_i, numNodes;
	int* 	nodes 	= bcs->GetSide( side, &numNodes );

	for( node_i = 0; node_i < numNodes; node_i++ ) {
		vals[nodes[node_i]][dof] = bcFunc( mesh->verts[nodes[node_i]] );
	}
}

void Field::SetICConst( int dof, double val ) {
	for( int node_i = 0; node_i < mesh->nVertsTotal; node_i++ ) {
		vals[node_i][dof] = val;
	}
}

void Field::SetICFunc( int dof, BCFunc* icFunc ) {
	for( int node_i = 0; node_i < mesh->nVertsTotal; node_i++ ) {
		vals[node_i][dof] = icFunc( mesh->verts[node_i] );
	}
}

void Field::PeriodicUpdate() {
	if( mesh->periodic[0] ) {
		for( int node_i = 0; node_i < mesh->nVerts[1]; node_i++ ) {
			for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
				vals[node_i*mesh->nVerts[0] + mesh->nVerts[0] - 1][dof_i] = vals[node_i*mesh->nVerts[0]][dof_i];
			}
		}
	}
	if( mesh->periodic[1] ) {
		for( int node_i = 0; node_i < mesh->nVerts[0]; node_i++ ) {
			for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
				vals[mesh->nVerts[0]*(mesh->nVerts[1] - 1) + node_i][dof_i] = vals[node_i][dof_i];
			}
		}
	}
	if( mesh->periodic[0] && mesh->periodic[1] ) {
		for( int dof_i = 0; dof_i < nDofs; dof_i++ ) {
			vals[mesh->nVerts[0]*mesh->nVerts[1] - 1][dof_i] = vals[0][dof_i];
		}
	}
}

//#define HDF5_OLD

void Field::Save( int timeStep ) {
	char		filename[50];
	hid_t		file, attribData_id, attrib_id, group_id, fileSpace, fileData, memSpace;
	hsize_t		a_dims, start[2], count[2], size[2];
	int		nDims	= 2;
	int		attribData;
	int		nLinearEls[2];
	double		buf[nDofs];
	int		node_i, dof_i;

	sprintf( filename, "%s.%.5u.h5", name.c_str(), timeStep );
	file = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

	/* write the dimensionality */
	a_dims = 1;
	attribData = nDims;
	attribData_id = H5Screate_simple( 1, &a_dims, NULL );
#ifdef HDF5_OLD
	group_id  = H5Gopen( file, "/" );
	attrib_id = H5Acreate( group_id, "dimensions", H5T_STD_I32BE, attribData_id, H5P_DEFAULT );
#else
	group_id  = H5Gopen( file, "/", H5P_DEFAULT );
	attrib_id = H5Acreate( group_id, "dimensions", H5T_STD_I32BE, attribData_id, H5P_DEFAULT, H5P_DEFAULT );
#endif
	H5Awrite( attrib_id, H5T_NATIVE_INT, &attribData );
	H5Aclose( attrib_id );
	H5Gclose( group_id );
	H5Sclose( attribData_id );
	
	/* store the mesh resolution */
	a_dims = nDims;
	nLinearEls[0] = mesh->nVerts[0] - 1;
	nLinearEls[1] = mesh->nVerts[1] - 1;
	attribData_id = H5Screate_simple( 1, &a_dims, NULL );
#ifdef HDF5_OLD 
	group_id  = H5Gopen(file, "/" );
	attrib_id = H5Acreate(group_id, "mesh resolution", H5T_STD_I32BE, attribData_id, H5P_DEFAULT );
#else
	group_id  = H5Gopen(file, "/", H5P_DEFAULT );
	attrib_id = H5Acreate(group_id, "mesh resolution", H5T_STD_I32BE, attribData_id, H5P_DEFAULT, H5P_DEFAULT );
#endif
        H5Awrite( attrib_id, H5T_NATIVE_INT, nLinearEls );
	H5Aclose( attrib_id );
	H5Gclose( group_id );
	H5Sclose( attribData_id );

	size[0] = mesh->nVertsTotal;
	size[1] = nDofs;
	fileSpace = H5Screate_simple( 2, size, NULL );
#ifdef HDF5_OLD
	fileData  = H5Dcreate( file, "/data", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT );
#else
	fileData  = H5Dcreate( file, "/data", H5T_NATIVE_DOUBLE, fileSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
#endif

	count[0] = 1;
	count[1] = nDofs;
	memSpace = H5Screate_simple( 2, count, NULL );
	H5Sselect_all( memSpace );
	for( node_i = 0; node_i < mesh->nVertsTotal; node_i++ ) {
		for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
			buf[dof_i] = vals[node_i][dof_i];
		}
		start[0] = node_i;
		start[1] = 0;
		H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Dwrite( fileData, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, buf );
	}

	H5Dclose( fileData );
	H5Sclose( memSpace );
	H5Sclose( fileSpace );
	H5Fclose( file );
}

void Field::Read( int timeStep ) {
	hid_t		file, fileData, fileSpace, memSpace;
	hsize_t		count[2], start[2];
	int		node_i, dof_i;
	double*		buf		= new double[nDofs];
	char		filename[50];

	sprintf( filename, "%s.%.5u.h5", name.c_str(), timeStep );

	file = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
#ifdef HDF5_OLD
	fileData = H5Dopen( file, "/data" );
#else
	fileData = H5Dopen( file, "/data", H5P_DEFAULT );
#endif
	fileSpace = H5Dget_space( fileData );

	count[0] = 1;
	count[1] = nDofs;
	memSpace = H5Screate_simple( 2, count, NULL );
	start[1] = 0;

	for( node_i = 0; node_i < mesh->nVertsTotal; node_i++ ) {
		start[0] = node_i;
		H5Sselect_hyperslab( fileSpace, H5S_SELECT_SET, start, NULL, count, NULL );
		H5Sselect_all( memSpace );
		H5Dread( fileData, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, buf );
		for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
			vals[node_i][dof_i] = buf[dof_i];
		}
	}
	delete[] buf;

	H5Dclose( fileData );
	H5Sclose( memSpace );
	H5Sclose( fileSpace );
	H5Fclose( file );
}
