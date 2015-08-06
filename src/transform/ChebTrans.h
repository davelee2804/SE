class ChebTrans {
	public:
		ChebTrans( int n, Mesh* mesh, int dim );
		~ChebTrans();
		double	Tn( int m, double x );
		void 	Transform( double* in );
		double*	Tj;
		double* x_gl;
		double*	x_gl_global;
	private:
		double*	T;
		double*	Tinv;
		int	n;
};
