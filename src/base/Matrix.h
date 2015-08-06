class Matrix {
	public:
		Matrix( std::string _name, Mat _mat, Vector* _rowVector, Vector* _colVector, Vector* _rhsVector, Operator** _ops, int _nOps );
		~Matrix();
		std::string name;
		Vector* rowVector;
		Vector* colVector;
		Vector* rhsVector;
		Operator** ops;
		int nOps;
		Mat mat;
		void AssembleElement( int el_i, double* M );
		void Assemble();
		void Write( std::string filename, bool writeZeros );
		void InverseDiagonal( Mat Minv );
};
