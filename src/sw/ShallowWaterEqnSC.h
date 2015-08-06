class ShallowWaterEqnSC {
	public:
		ShallowWaterEqnSC( Field* _vel, Field* _eta, double _dt, ShallowWaterParams* _params, Field* _topo, Field* _velBar, int svvCutoff );
		~ShallowWaterEqnSC();
		Field* 			vel;
		Field* 			eta;
		double 			dt;
		ShallowWaterParams* 	params;
		Field*			topo;
		Field*			velBar;
		void			Solve( int order, Field* velFieldRhs );
		void			Assemble( int order );
	private:
		SchurComplement* 	sc;
		HelmholtzEqn*		hh;
		Mat			Kinv;
		Mat			G;
		Mat			D;
		Mat			C;
		Vec			b;
};
