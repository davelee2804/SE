/* specreal element mesh - takes as arguments the number of elements in each dimension, and the order of the elements
 * in each dimension and generates a regular cartesian mesh of elements which chebychev basis functions, and quadrature
 * points at the roots of the associanted legendre polynomials */
class Mesh {
	public:
		Mesh( std::string _name, int* _nEls, std::string elType, int _N, double* _min, double* _max, bool* _periodic );
		~Mesh();
		std::string name;
		Element* el;
		int dim;
		int nVertsTotal;
		int nElsTotal;
		int nVerts[2];
		int nEls[2];
		double dx[2];
		double min[2];
		double max[2];
		bool periodic[2];
		double** verts;
		double** gGNix;
		double** gGNixx;
		int* elNodes;
		int N;
		void GetEl( double* coord, int* el_i );
		int* ElNodes( int el_i );
		void IndexToTopo( int index, int* topo );
		void TopoToIndex( int* topo, int* index );
		void GlobalToLocal( double* gCoord, int* el_i, double* lCoord );
		void LocalToGlobal( double* lCoord, int el_i, double* gCoord );
		double DetJac( int el_i, int pt_i );
		double DetJacAtCoord( int el_i, double* coord );
		double** ShapeFuncDerivs( int el_i, int pt_i, double* detJac );
		double** ShapeFuncDerivsAtCoord( int el_i, double* coord, double* detJac );
		double** ShapeFuncSecondDerivs( int el_i, int pt_i, double* detJac );
		double** ShapeFuncSecondDerivsAtCoord( int el_i, double* coord, double* detJac );
		void Save();
};
