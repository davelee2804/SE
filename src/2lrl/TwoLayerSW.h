typedef struct {
	double	p;
	double	g;
	double	f;
	double	beta;
	double	nu;
	double	H1;
	double	H2;
	double	tau;
	double	k;
} TwoLayerParams;

class TwoLayerSW {
	public:
		TwoLayerSW( Field* _velTop, Field* _velBot, Field* _etaInt, Field* _presSurf, TwoLayerParams* _params, double _dt );
		~TwoLayerSW();
		bool		pnull;
		bool		enull;
		Field*		velTop;
		Field*		velBot;
		Field*		etaInt;
		Field*		presSurf;
		TwoLayerParams*	params;
		void Solve( Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev );
	private:
		double		dt;
		bool		ass1;
		bool		ass2;
		Mat		PS;
		Mat		PI;
		Mat		UI;
		Mat		VI;
		void		AssembleMatrices( double a );
		void		AssembleGFirstOrder( Vec* G );
		void		AssembleGSecondOrder( Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev, Vec* G );
		void		CalcVelBar( Field* vel1, Field* vel2, Field* eta, Field* velb );
		void		SolvePressureFirstOrder();
		void		SolvePressureSecondOrder( Field* velTopPrev, Field* velBarPrev, Field* etaIntPrev );
		void		SolveFirstOrder();
		void		SolveSecondOrder( Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev );
};
