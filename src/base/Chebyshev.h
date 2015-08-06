class Chebyshev : public Element {
	public:
		Chebyshev( int _N );
		virtual ~Chebyshev();
		virtual double*  ShapeFuncs( int pt );
		virtual double** ShapeFuncDerivs( int pt );
		virtual double*  ShapeFuncsAtCoord( double* x );
		virtual double** ShapeFuncDerivsAtCoord( double* x );
		double** ShapeFunc2ndDerivs( int pt );
	private:
		double*  	dCij;
		double*  	d2Cij;
		double** 	GNixx;
		Jacobi**	polys;
		double		Cj( int j, double x );
		double		dCjdx( int j, double x );
};
