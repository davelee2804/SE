class Legendre : public Element {
	public:
		Legendre( int _N );
		virtual ~Legendre();
		virtual double* 	ShapeFuncs( int pt );
		virtual double** 	ShapeFuncDerivs( int pt );
		virtual double* 	ShapeFuncsAtCoord( double* coord );
		virtual double** 	ShapeFuncDerivsAtCoord( double* coord );
		double** 		ShapeFuncSecondDerivs( int i );
		double**		ShapeFuncSecondDerivsAtCoord( double* x );
		Jacobi** 		ModalBasis();
		double*  		ModalToNodalTransformMatrix();
	private:
		bool IsRoot( double x, int* i );
		double* dCij;
		double* d2Cij;
		double** GNixx;
};
