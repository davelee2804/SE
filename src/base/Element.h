class Element {
	public:
		Element( int _N );
		virtual ~Element();
		int N;
		int nNodes;
		int nPoints;
		QuadPoint** quadPts;
		Jacobi* P;
		Jacobi** polys;
		double** xi;
		double* abcissa;
		double* weights;
		double* Ni;
		double** GNix;
		virtual double* ShapeFuncs( int pt ) = 0;
		virtual double** ShapeFuncDerivs( int pt ) = 0;
		virtual double* ShapeFuncsAtCoord( double* coord ) = 0;
		virtual double** ShapeFuncDerivsAtCoord( double* coord ) = 0;
};

class Lagrange: public Element {
	public:
		Lagrange( int _N );
		virtual ~Lagrange();
		virtual double* ShapeFuncs( int pt );
		virtual double** ShapeFuncDerivs( int pt );
		virtual double* ShapeFuncsAtCoord( double* coord );
		virtual double** ShapeFuncDerivsAtCoord( double* coord );
		double Cj( int j, double xj, double x );
};

class Trig : public Element {
	public:
		Trig( int _N );
		virtual ~Trig();
		virtual double* ShapeFuncs( int pt );
		virtual double** ShapeFuncDerivs( int pt );
		virtual double* ShapeFuncsAtCoord( double* coord );
		virtual double** ShapeFuncDerivsAtCoord( double* coord );
		double Cj( int j, double xj, double x );
		double dCjdx( int j, double xj, double x );
	private:
		double Ninv;
};

class Linear : public Element {
	public:
		Linear( int _N );
		virtual ~Linear();
		virtual double* ShapeFuncs( int pt );
		virtual double** ShapeFuncDerivs( int pt );
		virtual double* ShapeFuncsAtCoord( double* coord );
		virtual double** ShapeFuncDerivsAtCoord( double* coord );
};

class Quadratic : public Element {
	public:
		Quadratic( int _N );
		virtual ~Quadratic();
		virtual double* ShapeFuncs( int pt );
		virtual double** ShapeFuncDerivs( int pt );
		virtual double* ShapeFuncsAtCoord( double* coord );
		virtual double** ShapeFuncDerivsAtCoord( double* coord );
};

double RootFinder( double x0, Jacobi* P );
bool   AlreadyFound( int n, double* roots, double root );
void   SortEntries( double* vals, int n );
