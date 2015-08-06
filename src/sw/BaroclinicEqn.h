class BaroclinicEqn {
	public:
		BaroclinicEqn( Field* 			_v1Curr, 
			       Field* 			_v1Prev, 
			       Field* 			_v2Curr, 
			       Field* 			_v2Prev, 
			       Field* 			_etaCurr, 
			       Field* 			_etaPrev, 
			       Field* 			_psCurr, 
			       Field* 			_psPrev, 
			       ShallowWaterParams* 	_topParams, 
			       ShallowWaterParams* 	_botParams, 
			       Field* 			_topo, 
			       int 			_svvCutoff,
			       double 			_dt );
		~BaroclinicEqn();
		void Solve( int order );
		Field* v1Curr;				/* top layer current velocity */
		Field* v1Prev;				/* top layer previous velocity */
		Field* v2Curr;				/* bottom layer current velocity */
		Field* v2Prev;				/* bottom layer previous velocity */
		Field* etaCurr;				/* current interface preturbation field */
		Field* etaPrev;				/* previous interface perturbation field */
		Field* psCurr;				/* current surface pressure */
		Field* psPrev;				/* previous surface pressure */
		Field* topo;				/* topography field */
		Field* btVel;				/* barotropic velocity */
		Field* bcVel1;				/* top layer baroclinic velocity */
		Field* bcVel2;				/* bottom layer baroclinic velocity */
		ShallowWaterParams* topParams;		/* top layer parameters */
		ShallowWaterParams* botParams;		/* bottom layer params */
	private:
		int 	svvCutoff;
		double	dt;
		SchurComplement* scSolver;		/* schur complement solver (for the implicit barotropic system) */
		ShallowWaterEqn* swSolver;		/* shallow water solver (for the bottom layer) */
		void AssembleBarotropic();
		void BaroclinicPredictorSetup();
		void BaroclinicPredictor( Field* bcVel1P, Field* bcVel2P, Field* etaP );
		void SolveBarotropic( Field* btv, Field* ps, Field* husq );
		void SolveShallowWater( Field* btv, Field* husq, int order );
		void CalcKineticEnergyFlux( Field* eta, Field* velocity1, Field* velocity2, Field* husq );
		void CalcHeightFields( Field* eta );
		Mat bt_Kinv;				/* barotropic system mass+coriolis inverse matrix */
		Mat bt_G;				/* barotropic system gradient matrix */
		Mat bt_D;				/* barotropic system divergence matrix */
		Vec bt_b;				/* barotropic system velocity boundary conditions vector */
		Mat bc_M_vp;				/* baroclinic predictor step velocity matrix */
		Mat bc_M_ep;				/* baroclinic predictor step interface perturbation matrix */
		Vec bc_b_v;				/* baroclinic predictor step velocity boundary conditions vector */
		Field* btVelP;				/* predictor step barotropic velocity */
		Field* bcVel1P;				/* predictor step top layer baroclinic velocity */
		Field* bcVel2P;				/* predictor step bottom layer baroclinic velocity */
		Field* etaP;				/* predictor step interface perturbation field */
		Field* husqN;				/* kinetic energy flux at prevous time level */
		Field* husqP;				/* kinetic energy flux at predictor time level */
		Field* topHeight;			/* top layer height field */
		Field* botHeight;			/* bottom layer height field */
		Advector* v1Adv;			/* top layer nonlinear advection operator advector */
		Advector* v2Adv;			/* bottom layer nonlinear advection operator advector */
		bool swAssembled;
};
