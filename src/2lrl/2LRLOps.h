class BetaInvMatrix : public Operator {
        public:
                BetaInvMatrix( std::string _name, Field* _rowField, Field* _colField, double _constant, double _alpha, double _f0, double _beta );
                virtual ~BetaInvMatrix();
		double alpha;
                double f0;
                double beta;
                virtual void AssembleElement( int el_i, double* M );
};

class PhiMatrix : public Operator {
        public:
                PhiMatrix( std::string _name, Field* _rowField, Field* _colField, double _constant, double _alpha, double _dt, double _f0, double _beta );
                virtual ~PhiMatrix();
		double alpha;
		double dt;
                double f0;
                double beta;
                virtual void AssembleElement( int el_i, double* M );
};

class GLayerVector : public RHSOp {
	public:
		GLayerVector( std::string 	_name, 
			      Mesh* 		_mesh, 
			      double 		_constant, 
			      Field* 		_field, 	/* topography */
			      Field* 		_height1, 
			      Field* 		_height2, 
			      Field* 		_velocity1, 
			      Field* 		_velocity2, 
			      double		_H1,
			      double		_H2,
			      double		_nu,
			      double		_gamma,
			      double		_gPrime,
			      double		_tau0,
			      double		_kws );
		virtual 	~GLayerVector();
		Field* 		height1;
		Field* 		height2;
		Field* 		velocity1;
		Field* 		velocity2;
		double		H1;
		double		H2;
		double		nu;
		double		gamma;
		double		gPrime;
		double		tau0;
		double		kws;
		double** 	db;
		double**	dh1;
		double**	dh2;
		double**	du1;
		double**	du2;
		virtual void AssembleElement( int el_i, double* G );
};

class GLayerVector2 : public RHSOp {
	public:
		GLayerVector2( std::string 	_name, 
			       Mesh* 		_mesh, 
			       double 		_constant, 
			       Field* 		_field,		/* eta-int */
			       Field* 		_velocity1, 
			       Field* 		_velocity2, 
			       double		_H1,
			       double		_H2,
			       double		_nu,
			       double		_gPrime,
			       double		_tau0,
			       double		_kws );
		virtual 	~GLayerVector2();
		Field* 		velocity1;
		Field* 		velocity2;
		double		H1;
		double		H2;
		double		nu;
		double		gPrime;
		double		tau0;
		double		kws;
		double**	dei;
		double**	du1;
		double**	du2;
		virtual void AssembleElement( int el_i, double* G );
};

class DivTotHeightVelRHS : public RHSOp {
	public:
		DivTotHeightVelRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field, Field* _velocity, double _H, Field* _topo );
		virtual ~DivTotHeightVelRHS();
		virtual void AssembleElement( int el_i, double* rhs );
	private:
		Field* 	velocity;
		Field*	topo;
		double 	H;
};

class BetaRHS : public RHSOp {
	public:
		BetaRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field, double _f0, double _beta );
		virtual ~BetaRHS();
		virtual void AssembleElement( int el_i, double* rhs );
	private:
		double f0;
		double beta;
};

class BaroclinicPressureBC : public RHSOp {
	public:
		BaroclinicPressureBC( std::string _name, Mesh* _mesh, double _constant, Field* _field, 
				      Params* _p, Field* _velTop, Field* _velBot, Field* _velBar, Field* _height1, Field* _height2 );
		virtual ~BaroclinicPressureBC();
		virtual void AssembleElement( int el_i, double* rhs );
	private:
		Params*		p;
		Field*		velTop;
		Field*		velBot;
		Field*		velBar;
		Field* 		height1;
		Field*		height2;
		double** 	dh;
		double**	du1;
		double**	du2;
		void 		AssembleBottom( int el_i, double* rhs );
		void 		AssembleTop( int el_i, double* rhs );
		void 		AssembleLeft( int el_i, double* rhs );
		void 		AssembleRight( int el_i, double* rhs );
};

class BarotropicPressureBC : public RHSOp {
	public:
		BarotropicPressureBC( std::string _name, Mesh* _mesh, double _constant, Field* _field );
		virtual ~BarotropicPressureBC();
		virtual void AssembleElement( int el_i, double* rhs );
	private:
		void AssembleBottom( int el_i, double* rhs );
		void AssembleTop( int el_i, double* rhs );
		void AssembleLeft( int el_i, double* rhs );
		void AssembleRight( int el_i, double* rhs );
};

class WindStress2RHS : public RHSOp {
	public:
		WindStress2RHS( std::string _name, Mesh* _mesh, double _constant, Field* _field, double _k, double _H );
		virtual ~WindStress2RHS();
		double k;
		double H;
		virtual void AssembleElement( int el_i, double* rhs );
};

class PhiRHS : public RHSOp {
	public:
		PhiRHS( std::string _name, Mesh* _mesh, double _constant, Field* field,
			double _sign, double _dt, double _H, double _tau0, double _k, double _nu, double _p, double _f0, double _beta,
			Field* _pres, Field* _eta, Field* _fieldSL );
		virtual ~PhiRHS();
		double		sign;
		double		dt;
		double 		H;
		double 		tau0;
		double 		k;
		double 		nu;
		double		p;
		double		f0;
		double		beta;
		Field* 		pres;
		Field* 		eta;
		Field* 		fieldSL;
		double** 	d2u;
		double**	dp;
		virtual void AssembleElement( int el_i, double* rhs );
};
