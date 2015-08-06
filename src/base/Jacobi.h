class Jacobi {
	public:
		Jacobi( int _n, double _a );
		~Jacobi();
		int n;
		double a;
		double* c;
		void Gen( Jacobi* m1, Jacobi* m2 );
		double Eval( double x );
		double EvalTotal( double x );
		double EvalDeriv( double x );
		double Eval2ndDeriv( double x );
};
