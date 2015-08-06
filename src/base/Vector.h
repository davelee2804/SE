class Vector {
	public:
		Vector( std::string _name, Field* _field, Vec _vec, RHSOp** _ops, int _nOps );
		~Vector();
		int vecSize;
		std::string name;
		Field* field;
		RHSOp** ops;
		int nOps;
		Vec vec;
		void UpdateField();
		void UpdateVec();
		void AssembleElement( int el_i, double* rhs );
		void Assemble();
		void CopyTo( Vec v );
		void CopyFrom( Vec v );
		void Write( std::string filename );
		double Norm();
};
