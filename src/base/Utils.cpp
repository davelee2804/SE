#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Legendre.h"
#include "Mesh.h"
#include "BCs.h"
#include "Field.h"
#include "Utils.h"

using namespace std;

/* NOTE: nodes are not evenly spaced for Chebychev polynomials, so this estimate is a little crude */
double CalcTimeStep( Field* velocity, double courantNo ) {
	double dx = velocity->mesh->dx[0];
	double dy = velocity->mesh->dx[1];
	double minSep = ( dx < dy ) ? dx : dy;
	double maxVel;
	int vert_i;
	double velSq, maxVelSq = 0.0;

	for( vert_i = 0; vert_i < velocity->mesh->nVertsTotal; vert_i++ ) {
		velSq = velocity->vals[vert_i][0]*velocity->vals[vert_i][0] + velocity->vals[vert_i][1]*velocity->vals[vert_i][1];
		if( velSq > maxVelSq ) {
			maxVelSq = velSq;
		}
	}
	maxVel = sqrt( maxVelSq );

	return courantNo*minSep/maxVel;
}

/* NOTE: nodes are not evenly spaced for Chebychev polynomials, so this estimate is a little crude */
double CalcViscosity( Mesh* mesh ) {
	double maxSep;
	
	//maxSep = (mesh->dx[0] > mesh->dx[1]) ? mesh->dx[0]/mesh->el->N : mesh->dx[1]/mesh->el->N;
	maxSep = (mesh->max[0] - mesh->min[0])/(mesh->nEls[0]*mesh->el->N);
	//maxSep = (mesh->max[0] - mesh->min[0])/(mesh->nEls[0]);
	//return maxSep*maxSep/(4.0*M_PI*M_PI);
	return maxSep*maxSep/(M_PI*M_PI);
}

void WriteXDMFHeader( int timeStep ) {
	ofstream 	file;
	char 		filename[20];
	char		ts[6] = "00000";

	cout << "writing fields to file at time step: " << timeStep << endl;

	sprintf( ts, "%.5u", timeStep );
	sprintf( filename, "XDMF.%s.xmf", ts );
	file.open( filename );

	file << "<?xml version=\"1.0\" ?>\n";
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
	file << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
	file << "\n";
	file << "<Domain>\n";
	file << "\n";

	file.close();
}

void WriteXDMFFooter( int timeStep ) {
	ofstream 	file;
	char 		filename[20];
	char		ts[6] = "00000";

	sprintf( ts, "%.5u", timeStep );
	sprintf( filename, "XDMF.%s.xmf", ts );
	file.open( filename, ios::app );

	file << "</Domain>\n";
	file << "\n";
	file << "</Xdmf>\n";
	file << "\n";

	file.close();
}

void WriteXDMF( Field** fields, int nFields, int timeStep, double time, double dt ) {
	ofstream 	file;
	char		filename[20];
	int 		field_i;
	Mesh*		mesh = fields[0]->mesh;
	char		varType[40];
	Field*		field;
	char		ts[6] = "00000";
	int		nLinearEls;

	sprintf( ts, "%.5u", timeStep );

	/* write the timestepping info to file */
	sprintf( filename, "time.%s.ts", ts );
	file.open( filename );
	file << "time step:\t" << timeStep << "\ttime:\t" << time << "\tdt:\t" << dt << "\n";
	file.close();

	for( field_i = 0; field_i < nFields; field_i++ ) {
		fields[field_i]->Save( timeStep );
	}

	sprintf( varType, "NumberType=\"Float\" Precision=\"8\"" );

	sprintf( filename, "XDMF.%s.xmf", ts );
	file.open( filename, ios::app );

	/* mesh */
	nLinearEls = (mesh->nVerts[0] - 1)*(mesh->nVerts[1] - 1);
        file << "   <Grid Name=\"FEM_Grid_" << mesh->name << "\">\n\n";
        file << "      <Time Value=\"" << time << "\" />\n\n";
        file << "         <Topology Type=\"Quadrilateral\" NumberOfElements=\"" << nLinearEls << "\"> \n";
        file << "            <DataItem Format=\"HDF\" DataType=\"Int\"  Dimensions=\"" << nLinearEls << " 4\">" << mesh->name << ".h5:/connectivity</DataItem>\n";
	file << "         </Topology>\n\n";
	file << "         <Geometry Type=\"XYZ\">\n";
	file << "            <DataItem ItemType=\"Function\"  Dimensions=\"" << mesh->nVertsTotal << " 3\" Function=\"JOIN($0, $1, 0*$1)\">\n";
	file << "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << mesh->nVertsTotal << " 1\" Name=\"XCoords\">\n";
	file << "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 " << mesh->nVertsTotal << " 1 </DataItem>\n";
	file << "                  <DataItem Format=\"HDF\" " << varType << " Dimensions=\"" << mesh->nVertsTotal << " 2\">" << mesh->name << ".h5:/vertices</DataItem>\n";
	file << "               </DataItem>\n";
	file << "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << mesh->nVertsTotal << " 1\" Name=\"YCoords\">\n";
	file << "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 1 1 1 " << mesh->nVertsTotal << " 1 </DataItem>\n";
	file << "                  <DataItem Format=\"HDF\" " << varType << " Dimensions=\"" << mesh->nVertsTotal << " 2\">" << mesh->name << ".h5:/vertices</DataItem>\n";
	file << "               </DataItem>\n";
	file << "            </DataItem>\n";
	file << "         </Geometry>\n\n";

	/* fields */
	for( field_i = 0; field_i < nFields; field_i++ ) {
		field = fields[field_i];
		if ( field->nDofs == 1 ) {
			file << "         <Attribute Type=\"Scalar\" Center=\"Node\" Name=\"" << field->name << "\">\n";
			file << "            <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << mesh->nVertsTotal << " 1\" >\n";
			file << "               <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 " << mesh->nVertsTotal << " 1 </DataItem>\n";
			file << "               <DataItem Format=\"HDF\" " << varType << " Dimensions=\"" << mesh->nVertsTotal << " 1\">" << field->name << "." << ts << ".h5:/data</DataItem>\n";
			file << "            </DataItem>\n";
			file << "         </Attribute>\n\n";
		} 
		else if( field->nDofs == 2 ) {
			file << "         <Attribute Type=\"Vector\" Center=\"Node\" Name=\"" << field->name << "\">\n";
			file << "            <DataItem ItemType=\"Function\"  Dimensions=\"" << mesh->nVertsTotal << " 3\" Function=\"JOIN($0, $1, 0*$1)\">\n";
			file << "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << mesh->nVertsTotal << " 1\" Name=\"XValue\">\n";
			file << "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 " << mesh->nVertsTotal << " 1 </DataItem>\n";
			file << "                  <DataItem Format=\"HDF\" " << varType << " Dimensions=\"" << mesh->nVertsTotal << " 2\">" << field->name << "." << ts << ".h5:/data</DataItem>\n";
			file << "               </DataItem>\n";
			file << "               <DataItem ItemType=\"HyperSlab\" Dimensions=\"" << mesh->nVertsTotal << " 1\" Name=\"YValue\">\n";
			file << "                  <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 1 1 1 " << mesh->nVertsTotal << " 1 </DataItem>\n";
			file << "                  <DataItem Format=\"HDF\" " << varType << " Dimensions=\"" << mesh->nVertsTotal << " 2\">" << field->name << "." << ts << ".h5:/data</DataItem>\n";
			file << "               </DataItem>\n";
			file << "            </DataItem>\n";
			file << "         </Attribute>\n\n";
		}
		else {
			cerr << "number of dofs cannot be written\n";
		}
	}
	file << "   </Grid>\n\n";

	file.close();
}

void WriteXDMFTemporal( int nTimeSteps, int dumpEvery ) {
	ofstream	file;
	char		filename[40];
	int		step_i;

	sprintf( filename, "XDMF.temporalAll.xmf" );
	file.open( filename );
	/* header info */
        file << "<?xml version=\"1.0\" ?>\n";
        file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        file << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
        file << "\n";
        file << "<Domain>\n";
        file << "\n";
	/* fields info */
	file << "   <xi:include href=\"XDMF.FilesField.xdmf\" xpointer=\"xpointer(//Xdmf/Grid)\"/>\n\n";
	/* footer info */
	file << "</Domain>\n";
        file << "\n";
        file << "</Xdmf>\n";
        file << "\n";
	file.close();

	sprintf( filename, "XDMF.temporalFields.xmf" );
	file.open( filename );
	/* header info */
        file << "<?xml version=\"1.0\" ?>\n";
        file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        file << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
        file << "\n";
        file << "<Domain>\n";
        file << "\n";
	/* fields info */
	file << "   <xi:include href=\"XDMF.FilesField.xdmf\" xpointer=\"xpointer(//Xdmf/Grid)\"/>\n\n";
	/* footer info */
	file << "</Domain>\n";
        file << "\n";
        file << "</Xdmf>\n";
        file << "\n";
	file.close();

	sprintf( filename, "XDMF.FilesField.xdmf" );
	file.open( filename );
	file << "<?xml version=\"1.0\" ?>\n";
	file << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n";
	file << "<Grid GridType=\"Collection\" CollectionType=\"Temporal\" Name=\"FEM_Mesh_Fields\">\n";

	for( step_i = 0; step_i <= nTimeSteps; step_i += dumpEvery ) {
		sprintf( filename, "XDMF.%.5u.xmf", step_i );
		file << "    <xi:include href=\"" << filename << "\" xpointer=\"xpointer(//Xdmf/Domain/Grid[1])\"/>\n";
	}

	file << "</Grid>\n";
        file << "</Xdmf>\n";

	file.close();
}

#define EPS 1.0e-6

double FieldError( Field* numeric, Field* analytic, int dof, bool normalise ) {
	double detJac, weight, *coord, *n, a[2], elErrSq, errSq = 0.0, gCoord[2], normSq = 0.0, elNormSq;
	int *elNodes;

	for( int el_i = 0; el_i < numeric->mesh->nElsTotal; el_i++ ) {
		elErrSq  = 0.0;
		elNormSq = 0.0;
		elNodes = numeric->mesh->ElNodes( el_i );
		for( int pt_i = 0; pt_i < numeric->mesh->el->nPoints; pt_i++ ) {
			coord  = numeric->mesh->el->quadPts[pt_i]->coord;
			weight = numeric->mesh->el->quadPts[pt_i]->weight;
			detJac = numeric->mesh->DetJac( el_i, pt_i );
			numeric->mesh->LocalToGlobal( coord, el_i, gCoord );
			n = numeric->vals[elNodes[pt_i]];
			analytic->InterpGlobal( gCoord, a );
			elErrSq  += detJac*weight*( n[dof] - a[dof] )*( n[dof] - a[dof] );
			elNormSq += detJac*weight*a[dof]*a[dof];
		}
		errSq  += elErrSq;
		normSq += elNormSq;
	}
	return (normalise) ? sqrt( errSq/normSq ) : sqrt( errSq );
}

void CreateVector( Field* field, Vec* v ) {
	int vSize;

	vSize = 0;
	for( int dof_i = 0; dof_i < field->nDofs; dof_i++ ) {
		vSize += field->mesh->nVertsTotal - field->bcs->size[dof_i];
	}

	VecCreate( MPI_COMM_WORLD, v );
	VecSetSizes( *v, vSize, PETSC_DETERMINE );
	VecSetType( *v, VECSEQ );
	VecSetOption( *v, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
	VecZeroEntries( *v );
}

void CreateMatrix( Field* rowField, Field* colField, Mat* A ) {
	int rowSize, colSize, allocSize;

	rowSize = 0;
	for( int dof_i = 0; dof_i < rowField->nDofs; dof_i++ ) {
		rowSize += rowField->mesh->nVertsTotal - rowField->bcs->size[dof_i];
	}
	colSize = 0;
	for( int dof_i = 0; dof_i < colField->nDofs; dof_i++ ) {
		colSize += colField->mesh->nVertsTotal - colField->bcs->size[dof_i];
	}

	allocSize = ( rowField->mesh->el->nNodes*rowField->nDofs > colField->mesh->el->nNodes*colField->nDofs ) ?
		      rowField->mesh->el->nNodes*rowField->nDofs : colField->mesh->el->nNodes*colField->nDofs;

	MatCreate( MPI_COMM_WORLD, A );
	MatSetSizes( *A, rowSize, colSize, PETSC_DETERMINE, PETSC_DETERMINE );
	MatSetType( *A, MATSEQAIJ );
	MatSeqAIJSetPreallocation( *A, 4*allocSize, PETSC_NULL );
	MatZeroEntries( *A );
}

/* NOTE: only for scalar fields at the moment */
/*void CreateNodalToModalMatrix( Field* field, Mat* N2M ) {
	Legendre*	l		= (Legendre*)field->mesh->el;
	double*		B		= NULL;
	double*		Binv		= new double[field->mesh->el->nNodes*field->mesh->el->nNodes];
	int*		map		= field->bcs->fieldToVecMap;
	int*		gNodes;
	int*		nodes		= NULL;

	if( field->nDofs > 1 ) {
		cout << "ERROR nodal->modal transform func only for scalar fields at present.\n";
		exit(0);
	}

	CreateMatrix( field, field, N2M );

	B = l->ModalToNodalTransformMatrix();
	Inv( B, Binv, field->mesh->el->nNodes );

	for( int el_i = 0; el_i < field->mesh->nElsTotal; el_i++ ) {
                gNodes = field->mesh->ElNodes( el_i );
                for( int node_i = 0; node_i < field->mesh->el->nNodes; node_i++ ) {
			nodes[node_i] = map[gNodes[node_i]];
                }
		MatSetValues( *N2M, field->mesh->el->nNodes, nodes, field->mesh->el->nNodes, nodes, Binv, ADD_VALUES );
	}

	delete[] B;
	delete[] Binv;
}*/
