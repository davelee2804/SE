class OIFS {
	public:
		OIFS( Operator* _op, Field* _phi0, Field* _phi1, Field* _psi0, Field* _psi1, double _dt, int _nCycles );
		~OIFS();
		void Solve();
		Field*		phiTilde;
	private:
		Operator* 	op;
		Field*		phi0;
		Field*		phi1;
		Field*		psi0;
		Field*		psi1;
		double		dt;
		int		nCycles;
		Mat		M;
		Mat		Minv;
		Mat		A;
		Vec		x;
		Vector*		vector;
		Matrix*		matrix;
		void		SolveFirstOrder();
		void		SolveSecondOrder();
		void		Extrapolate( Field* field, double a );
};
		
void InitMat( Mat* A, int rSize, int cSize, int alloc );
void InitVec( Vec* v, int size );
void FieldToVec( Field* field, Vec vec );
void VecToField( Vec vec, Field* field );
void FieldCopy( Field* from, Field* to );
