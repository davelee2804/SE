class SW2L {
	public:
		SW2L( Field* _velTop, Field* _velBot, Field* _etaInt, 
		      Field* _velTopPrev, Field* _velBotPrev, Field* _etaIntPrev, 
		      Field* _presSurf, Field* _bTopog, Params* _tp, Params* _bp, double _dt );
		~SW2L();
		bool		unull;
		bool		vnull;
		bool		enull;
		void		Assemble( int order );
		void		Solve( int order );
		Field*		velTop;
		Field*		velBot;
		Field*		etaInt;
		Params*		tp;
		Params*		bp;
	private:
		Field*		velTopPrev;
		Field*		velBotPrev;
		Field*		etaIntPrev;
		Field*		presSurf;
		Field*		bTopog;
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
