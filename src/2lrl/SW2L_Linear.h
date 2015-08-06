class SW2L_Linear {
	public:
		SW2L_Linear( Field* _velTop, Field* _velBot, Field* _etaInt, 
		      Field* _velTopPrev, Field* _velBotPrev, Field* _etaIntPrev, 
		      Field* _presSurf, Field* _bTopog, Params* _tp, Params* _bp, double _dt );
		~SW2L_Linear();
		bool		unull;
		bool		vnull;
		bool		enull;
		void		Assemble( int order );
		void		Solve( int order );
	private:
		Field*		velTop;
		Field*		velBot;
		Field*		etaInt;
		Field*		velTopPrev;
		Field*		velBotPrev;
		Field*		etaIntPrev;
		Field*		presSurf;
		Field*		bTopog;
		Params*		tp;
		Params*		bp;
		double		dt;
		Field*		velTopTemp;
		Field*		velBotTemp;
		Field*		etaIntTemp;
		Mat		A;
		Vec		b;
		void		FirstOrder();
		void		SecondOrder();
		void		SolveLinSys( Vec f, Vec x );
};

class GLayerVector_Linear : public RHSOp {
	public:
		GLayerVector_Linear( std::string 	_name, 
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
			      double		_kws,
			      Field*		_eta );
		virtual 	~GLayerVector_Linear();
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
		Field*		eta;
		double** 	db;
		double**	dh1;
		double**	dh2;
		double**	du1;
		double**	du2;
		double**	de;
		virtual void AssembleElement( int el_i, double* G );
};

void AssembleGFirstOrder_Linear( Params* tp, Params* bp, BarotropicEqn* bt, Field* velTop, Field* velBot, Field* etaInt, Vec* G );
void AssembleGSecondOrder_Linear( Params* tp, Params* bp, BarotropicEqn* bt, Field* velTop, Field* velBot, Field* etaInt, Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev, Vec* G );
void CalcVelBar_Linear( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* velBar );
