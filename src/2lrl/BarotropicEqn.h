class BarotropicEqn {
	public:
		BarotropicEqn( Field* _velocity, Field* _pressure, double _dt, Params* _params, Field* _topo );
		~BarotropicEqn();
		Field* 		velocity;
		Field* 		pressure;
		double 		dt;
		Params*		params;
		Field*		topo;
		bool		pnull;
		void 		Setup( int order );
		void 		Solve( Field* velPrev, Vec G );
		void		AssembleGFirstOrder( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* etaInt, Vec* G );
		void		AssembleGSecondOrder( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* etaInt, Field* velTopPrev, Field* velBotPrev, Field* etaIntPrev, Vec* G );
		void		CalcVelBar( Params* tp, Params* bp, Field* velTop, Field* velBot, Field* etaInt, Field* velBar );
		void		CalcHeights( Params* tp, Params* bp, Field* eta, Field** height1, Field** height2 );
		/*void		AssemblePressureBC( Params* tp, Params* bp, Params* bt, 
						    Field* velTop, Field* velTopPrev, 
						    Field* velBot, Field* velBotPrev, 
						    Field* etaInt, Field* etaIntPrev, 
						    Field* velBar, Field* velBarPrev, Vec* n );
		*/
		void		AssemblePressureBC( Field* velBar, Vec* n );
	private:
		Mat		M_u;
		Mat		M_p;
		Vec		b_u;
		Field*		totalHeight;
};
