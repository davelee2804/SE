#include <stdlib.h>
#include <iostream>
#include <string>
#include "QuadPoint.h"
#include "LinAlg.h"
#include "Jacobi.h"
#include "Element.h"
#include "Mesh.h"
#include "BCs.h"

using namespace std;
using std::string;

BCs::BCs( bool _bcBottom, bool _bcTop, bool _bcLeft, bool _bcRight, Mesh* _mesh, int _vecOffset ) {
	int vecNode, fieldNode, topo[2], node_i;

	mesh = _mesh;
	vecOffset = _vecOffset;

	bcBottom[0] = _bcBottom;
	bcTop[0] = _bcTop;
	bcLeft[0] = _bcLeft;
	bcRight[0] = _bcRight;

	bottom = top = left = right = NULL;

	size[0] = size[1] = 0;
	
	if( ( bcBottom[0] || bcTop[0] ) && mesh->periodic[1] ) {
		cerr << "attempting to assign BOTH periodic and non-periodic bcs in the vertical direction\n";
		exit(1);
	}
	if( ( bcLeft[0] || bcRight[0] ) && mesh->periodic[0] ) {
		cerr << "attempting to assign BOTH periodic and non-periodic bcs in the horizontal direction\n";
		exit(1);
	}

	if( mesh->periodic[0] ) {
		bcRight[0] = true;
	}
	if( mesh->periodic[1] ) {
		bcTop[0] = true;
	}
	if( bcBottom[0] ) {
		size[0] += mesh->nVerts[0];
		bottom = new int[mesh->nVerts[0]];
		for( node_i = 0; node_i < mesh->nVerts[0]; node_i++ ) {
			bottom[node_i] = node_i;
		}
	}
	if( bcTop[0] || mesh->periodic[1] ) {
		size[0] += mesh->nVerts[0];
		top = new int[mesh->nVerts[0]];
		for( node_i = 0; node_i < mesh->nVerts[0]; node_i++ ) {
			top[node_i] = mesh->nVerts[0]*(mesh->nVerts[1] - 1) + node_i;
		}
	}
	if( bcLeft[0] ) {
		size[0] += mesh->nVerts[1];
		left = new int[mesh->nVerts[1]];
		for( node_i = 0; node_i < mesh->nVerts[1]; node_i++ ) {
			left[node_i] = mesh->nVerts[0]*node_i;
		}
	}
	if( bcRight[0] || mesh->periodic[0] ) {
		size[0] += mesh->nVerts[1];
		right = new int[mesh->nVerts[1]];
		for( node_i = 0; node_i < mesh->nVerts[1]; node_i++ ) {
			right[node_i] = mesh->nVerts[0]*node_i + mesh->nVerts[0] - 1;
		}
	}
	if( bcBottom[0] && bcLeft[0] ) {
		size[0]--;
	}
	if( bcBottom[0] && bcRight[0] ) {
		size[0]--;
	}
	if( bcTop[0] && bcLeft[0] ) {
		size[0]--;
	}
	if( bcTop[0] && bcRight[0] ) {
		size[0]--;
	}

	vecToFieldMap = new int[mesh->nVertsTotal - size[0]];
	fieldToVecMap = new int[mesh->nVertsTotal];

	/* map petsc vector node indices to field node indices. if a field node isn't a boundary node, then 
	 * it is represented in the vector. */
	vecNode = 0;
	if( ( !mesh->periodic[0] && !mesh->periodic[1] ) || ( mesh->periodic[0] && mesh->periodic[1] ) ) {
		for( fieldNode = 0; fieldNode < mesh->nVertsTotal; fieldNode++ ) {
			mesh->IndexToTopo( fieldNode, topo );

			if( !IsBCNode( fieldNode ) ) {
				vecToFieldMap[vecNode] = fieldNode;
				fieldToVecMap[fieldNode] = vecNode + vecOffset;
				vecNode++;
			}
			else if( mesh->periodic[0] && mesh->periodic[1] && fieldNode == mesh->nVertsTotal - 1 ) {
				fieldToVecMap[fieldNode] = 0 + vecOffset;
			}
			else if( mesh->periodic[0] && topo[0] == mesh->nVerts[0] - 1 ) {
				fieldToVecMap[fieldNode] = topo[1]*(mesh->nVerts[0] - 1) + vecOffset;
			}
			else if( mesh->periodic[1] && topo[1] == mesh->nVerts[1] - 1 ) {
				fieldToVecMap[fieldNode] = topo[0] + vecOffset;
			}
			else {
				fieldToVecMap[fieldNode] = -1;
			}
		}
	}
	else if( mesh->periodic[0] && ( bcTop[0] && bcBottom[0] ) ) {
		for( fieldNode = 0; fieldNode < mesh->nVertsTotal; fieldNode++ ) {
			mesh->IndexToTopo( fieldNode, topo );

			if( !IsBCNode( fieldNode ) ) {
				vecToFieldMap[vecNode] = fieldNode;
				fieldToVecMap[fieldNode] = vecNode + vecOffset;
				vecNode++;
			}
			else if( topo[1] == 0 || topo[1] == mesh->nVerts[1] - 1 ) {
				fieldToVecMap[fieldNode] = -1;
			}
			else if( topo[0] == mesh->nVerts[0] - 1 ) { /* periodic in x */
				fieldToVecMap[fieldNode] = (topo[1] - 1)*(mesh->nVerts[0] - 1) + vecOffset;
			}
			else {
				cerr << "bc error: " << fieldNode << endl;
				exit(1);
			}
		}
	}
	else if( mesh->periodic[0] && bcTop[0] ) {
		for( fieldNode = 0; fieldNode < mesh->nVertsTotal; fieldNode++ ) {
			mesh->IndexToTopo( fieldNode, topo );

			if( !IsBCNode( fieldNode ) ) {
				vecToFieldMap[vecNode]   = fieldNode;
				fieldToVecMap[fieldNode] = vecNode + vecOffset;
				vecNode++;
			}
			else if( topo[1] == mesh->nVerts[1] - 1 ) { /* top bc */
				fieldToVecMap[fieldNode] = -1;
			}
			else if( topo[0] == mesh->nVerts[0] - 1 ) { /* periodic in x */
				fieldToVecMap[fieldNode] = topo[1]*(mesh->nVerts[0] - 1) + vecOffset;
			}
			else {
				cerr << "bc error: " << fieldNode << endl;
				exit(1);
			}
		}
	}
	else if( mesh->periodic[0] ) {
		for( fieldNode = 0; fieldNode < mesh->nVertsTotal; fieldNode++ ) {
			mesh->IndexToTopo( fieldNode, topo );
			if( !IsBCNode( fieldNode ) ) {
				vecToFieldMap[vecNode] = fieldNode;
				fieldToVecMap[fieldNode] = vecNode + vecOffset;
				vecNode++;
			}
			else if( topo[0] == mesh->nVerts[0] - 1 ) {
				fieldToVecMap[fieldNode] = topo[1]*(mesh->nVerts[0] - 1) + vecOffset;
			}
			else {
				cerr << "bc error: " << fieldNode << endl;
				exit(1);
			}
		}
	}
}

BCs::BCs( bool* _bcBottom, bool* _bcTop, bool* _bcLeft, bool* _bcRight, Mesh* _mesh, int _vecOffset ) {
	int dof_i, node_i, vecNode, fieldNode, fieldNodeAndDof, topo[2], sizeTotal;

	mesh = _mesh;
	vecOffset = _vecOffset;

	for( dof_i = 0; dof_i < 2; dof_i++ ) {
		bcBottom[dof_i] = _bcBottom[dof_i];
		bcTop[dof_i] = _bcTop[dof_i];
		bcLeft[dof_i] = _bcLeft[dof_i];
		bcRight[dof_i] = _bcRight[dof_i];
	}

	size[0] = size[1] = 0;

	if( ( bcBottom[0] || bcBottom[1] || bcTop[0] || bcTop[1] ) && mesh->periodic[1] ) {
		cerr << "attempting to assign BOTH periodic and non-periodic bcs in the vertical direction\n";
		exit(1);
	}
	if( ( bcLeft[0] || bcLeft[1] || bcRight[0] || bcRight[1] ) && mesh->periodic[0] ) {
		cerr << "attempting to assign BOTH periodic and non-periodic bcs in the horizontal direction\n";
		exit(1);
	}

	if( mesh->periodic[0] ) {
		bcRight[0] = true;
		bcRight[1] = true;
	}
	if( mesh->periodic[1] ) {
		bcTop[0] = true;
		bcTop[1] = true;
	}
	for( dof_i = 0; dof_i < 2; dof_i++ ) {
		if( bcBottom[dof_i] ) {
			size[dof_i] += mesh->nVerts[0];
		}
		if( bcTop[dof_i] ) {
			size[dof_i] += mesh->nVerts[0];
		}
		if( bcLeft[dof_i] ) {
			size[dof_i] += mesh->nVerts[1];
		}
		if( bcRight[dof_i] ) {
			size[dof_i] += mesh->nVerts[1];
		}
	}
	for( dof_i = 0; dof_i < 2; dof_i++ ) {
		if( bcBottom[dof_i] && bcLeft[dof_i] ) {
			size[dof_i]--;
		}
		if( bcBottom[dof_i] && bcRight[dof_i] ) {
			size[dof_i]--;
		}
		if( bcTop[dof_i] && bcLeft[dof_i] ) {
			size[dof_i]--;
		}
		if( bcTop[dof_i] && bcRight[dof_i] ) {
			size[dof_i]--;
		}
	}

	if( bcBottom[0] || bcBottom[1] ) {
		bottom = new int[mesh->nVerts[0]];
		for( node_i = 0; node_i < mesh->nVerts[0]; node_i++ ) {
			bottom[node_i] = node_i;
		}
	}
	if( bcTop[0] || bcTop[1] ) {
		top = new int[mesh->nVerts[0]];
		for( node_i = 0; node_i < mesh->nVerts[0]; node_i++ ) {
			top[node_i] = mesh->nVerts[0]*(mesh->nVerts[1] - 1) + node_i;
		}
	}
	if( bcLeft[0] || bcLeft[1] ) {
		left = new int[mesh->nVerts[1]];
		for( node_i = 0; node_i < mesh->nVerts[1]; node_i++ ) {
			left[node_i] = mesh->nVerts[0]*node_i;
		}
	}
	if( bcRight[0] || bcRight[1] ) {
		right = new int[mesh->nVerts[1]];
		for( node_i = 0; node_i < mesh->nVerts[1]; node_i++ ) {
			right[node_i] = mesh->nVerts[0]*node_i + mesh->nVerts[0] - 1;
		}
	}

	vecToFieldMap = new int[mesh->nVertsTotal*2 - size[0] - size[1]];
	fieldToVecMap = new int[mesh->nVertsTotal*2];

	vecNode = 0;
	sizeTotal = 0;
	/* no or full periodicity */
	if( ( !mesh->periodic[0] && !mesh->periodic[1] ) || ( mesh->periodic[0] && mesh->periodic[1] ) ) {
		for( dof_i = 0; dof_i < 2; dof_i++ ) {
			for( fieldNode = 0; fieldNode < mesh->nVertsTotal; fieldNode++ ) {
				fieldNodeAndDof = dof_i*mesh->nVertsTotal + fieldNode;
				mesh->IndexToTopo( fieldNode, topo );

				if( !IsBCNodeAndDof( fieldNode, dof_i ) ) {
					fieldToVecMap[fieldNodeAndDof] = vecNode + vecOffset;
					vecToFieldMap[vecNode] = dof_i*mesh->nVertsTotal + fieldNode;
					vecNode++;
				}
				else if( mesh->periodic[0] && mesh->periodic[1] && fieldNode == mesh->nVertsTotal - 1 ) {
					fieldToVecMap[fieldNodeAndDof] = dof_i*mesh->nVertsTotal - sizeTotal + vecOffset;
				}
				else if( mesh->periodic[0] && topo[0] == mesh->nVerts[0] - 1 ) {
					fieldToVecMap[fieldNodeAndDof] = dof_i*mesh->nVertsTotal - sizeTotal + topo[1]*(mesh->nVerts[0] - 1) + vecOffset;
				}
				else if( mesh->periodic[1] && topo[1] == mesh->nVerts[1] - 1 ) {
					fieldToVecMap[fieldNodeAndDof] = dof_i*mesh->nVertsTotal - sizeTotal + topo[0] + vecOffset;
				}
				else {
					fieldToVecMap[fieldNodeAndDof] = -1;
				}
			}
			sizeTotal += size[dof_i];
		}
	}
	/* periodic channel */
	else if( mesh->periodic[0] ) {
		for( dof_i = 0; dof_i < 2; dof_i++ ) {
			for( fieldNode = 0; fieldNode < mesh->nVertsTotal; fieldNode++ ) {
				fieldNodeAndDof = dof_i*mesh->nVertsTotal + fieldNode;
				mesh->IndexToTopo( fieldNode, topo );

				if( !IsBCNodeAndDof( fieldNode, dof_i ) ) {
					fieldToVecMap[fieldNodeAndDof] = vecNode + vecOffset;
					vecToFieldMap[vecNode] = dof_i*mesh->nVertsTotal + fieldNode;
					vecNode++;
				}
				else if( topo[1] == 0 && bcBottom[dof_i] ) {
					fieldToVecMap[fieldNodeAndDof] = -1;
				}
				else if( topo[1] == mesh->nVerts[1] - 1 && bcTop[dof_i] ) {
					fieldToVecMap[fieldNodeAndDof] = -1;
				}
				else if( topo[0] == mesh->nVerts[0] - 1 && bcBottom[dof_i] ) {
					fieldToVecMap[fieldNodeAndDof] = dof_i*mesh->nVertsTotal - sizeTotal + (topo[1] - 1)*(mesh->nVerts[0] - 1) + vecOffset;
				}
				else if( topo[0] == mesh->nVerts[0] - 1 ) {
					fieldToVecMap[fieldNodeAndDof] = dof_i*mesh->nVertsTotal - sizeTotal + topo[1]*(mesh->nVerts[0] - 1) + vecOffset;
				}
				else {
					fieldToVecMap[fieldNodeAndDof] = -1;
				}
			}
			sizeTotal += size[dof_i];
		}
	}
	else {
		cerr << "bc configuration not implemented\n";
		exit(1);
	}
}

BCs::~BCs() {
	delete[] vecToFieldMap;
	delete[] fieldToVecMap;

	if( bottom ) {
		delete[] bottom;
	}
	if( top ) {
		delete[] top;
	}
	if( left ) {
		delete[] left;
	}
	if( right ) {
		delete[] right;
	}
}

bool BCs::IsBCNode( int node ) {
	int topo[2];

	mesh->IndexToTopo( node, topo );
	if( bcBottom[0] && topo[1] == 0 ) {
		return true;
	}
	if( bcTop[0] && topo[1] == mesh->nVerts[1] - 1 ) {
		return true;
	}
	if( bcLeft[0] && topo[0] == 0 ) {
		return true;
	}
	if( bcRight[0] && topo[0] == mesh->nVerts[0] - 1 ) {
		return true;
	}

	return false;
}

bool BCs::IsBCNodeAndDof( int node, int dof ) {
	int topo[2];

	mesh->IndexToTopo( node, topo );
	if( topo[1] == 0 && bcBottom[dof] ) {
		return true;
	}
	if( topo[1] == mesh->nVerts[1] - 1 && bcTop[dof] ) {
		return true;
	}
	if( topo[0] == 0 && bcLeft[dof] ) {
		return true;
	}
	if( topo[0] == mesh->nVerts[0] - 1 && bcRight[dof] ) {
		return true;
	}

	return false;
}

int* BCs::GetSide( string side, int* nNodes ) {
	if( side == "bottom" ) {
		*nNodes = mesh->nVerts[0];
		return bottom;
	}
	if( side == "top" ) {
		*nNodes = mesh->nVerts[0];
		return top;
	}
	if( side == "left" ) {
		*nNodes = mesh->nVerts[1];
		return left;
	}
	if( side == "right" ) {
		*nNodes = mesh->nVerts[1];
		return right;
	}
	*nNodes = -1;
	return NULL;
}
