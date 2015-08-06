class Operator {
	public:
		Operator( std::string _name, Field* _rowField, Field* _colField, double _constant );
		virtual ~Operator();
		std::string name;
		Field* rowField;
		Field* colField;
		Field* field;
		double constant;
		virtual void AssembleElement( int el_i, double* M ) = 0;
};

class Laplacian : public Operator {
	public:
		Laplacian( std::string _name, Field* _rowField, Field* _colField, double _constant );
		virtual ~Laplacian();
		virtual void AssembleElement( int el_i, double* M );
};

class NonLinearLaplacian : public Operator {
	public:
		NonLinearLaplacian( std::string _name, Field* _rowField, Field* _colField, double _constant, Field* _field );
		virtual ~NonLinearLaplacian();
		Field* field;
		virtual void AssembleElement( int el_i, double* M );
};

class Biharmonic : public Operator {
        public:
                Biharmonic( std::string _name, Field* _rowField, Field* _colField, double _constant );
                virtual ~Biharmonic();
                virtual void AssembleElement( int el_i, double* M );
};

class StressTensor : public Operator {
	public:
		StressTensor( std::string _name, Field* _rowField, Field* _colField, double _constant );
		virtual ~StressTensor();
		virtual void AssembleElement( int el_i, double* M );
};

class MassMatrix : public Operator {
	public:
		MassMatrix( std::string _name, Field* _rowField, Field* _colField, double _constant );
		virtual ~MassMatrix();
		virtual void AssembleElement( int el_i, double* M );
};

class Gradient : public Operator {
	public:
		Gradient( std::string _name, Field* _rowField, Field* _colField, double _constant );
		virtual ~Gradient();
		virtual void AssembleElement( int el_i, double* M );
};

class Divergence : public Operator {
	public:
		Divergence( std::string _name, Field* _rowField, Field* _colField, double _constant );
		virtual ~Divergence();
		virtual void AssembleElement( int el_i, double* M );
};

class StokesPC : public Operator {
	public:
		StokesPC( std::string _name, Field* _rowField, Field* _colField, double _constant, Field* _velocity, Field* _pressure );
		virtual ~StokesPC();
		Field* velocity;
		Field* pressure;
		double* K;
		double* G;
		double* KinvG;
		Laplacian* laplacian;
		Gradient* gradient;
		virtual void AssembleElement( int el_i, double* M );
};

class BetaMatrix : public Operator {
	public:
		BetaMatrix( std::string _name, Field* _rowField, Field* _colField, double _constant, double _f0, double _beta );
		virtual ~BetaMatrix();
		double f0;
		double beta;
		virtual void AssembleElement( int el_i, double* M );
};

class ConvectionMatrix : public Operator {
        public:
                ConvectionMatrix( std::string _name, Field* _rowField, Field* _colField, double _constant, Field* _field );
                virtual ~ConvectionMatrix();
                virtual void AssembleElement( int el_i, double* M );
};
