class DivPhiVelMatrix : public Operator {
	public:
		DivPhiVelMatrix( std::string _name, Field* _rowField, Field* _colField, double _constant, Field* _field );
		virtual ~DivPhiVelMatrix();
		virtual void AssembleElement( int el_i, double* M );
};

class ConvMat1D : public Operator {
	public:
		ConvMat1D( std::string _name, Field* _rowField, Field* _colFIeld, double _constant, Field* _field );
		virtual ~ConvMat1D();
		virtual void AssembleElement( int el_i, double* M );
};

class WindStressRHS : public RHSOp {
	public:
		WindStressRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field, double _k, double _H );
		virtual ~WindStressRHS();
		double k;
		double H;
		virtual void AssembleElement( int el_i, double* rhs );
};

class DivHeightMinusTopo : public Operator {
	public:
		DivHeightMinusTopo( std::string _name, Field* _rowField, Field* _colField, double _constant, double _H, Field* _topo );
		virtual ~DivHeightMinusTopo();
		double H;
		Field* topo;
		double** db;
		virtual void AssembleElement( int el_i, double* rhs );
};

class DivVelBarMatrix : public Operator {
	public:
		DivVelBarMatrix( std::string _name, Field* _rowField, Field* _colField, double _constant, Field* _velBar );
		virtual ~DivVelBarMatrix();
		Field* velBar;
		virtual void AssembleElement( int el_i, double* M );
};

class DivVelBarRHS : public RHSOp {
	public:
		DivVelBarRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field, Field* _topo, double _H_i );
		virtual ~DivVelBarRHS();
		double H_i;
		Field* topo;
		virtual void AssembleElement( int el_i, double* rhs );
};
