class NavierStokesEqn {
	public:
		NavierStokesEqn( Field* velocity, Field* pressure, double _dt );
		~NavierStokesEqn();
		Field* velStep2;
		Field* divVelStep1;
		double dt;
		bool   firstOrderInit;
		void Solve( Field* velocity, Field* pressure, Field* prevVel, double nu, int svvCutoff );
	private:
		PoissonEqn* 	poisson;
		HelmholtzEqn* 	helmholtz;
		bool		secondOrderDone;
		void InitMats( int vSize, int pSize, int nVelNodes, int nPresNodes );
		void InitVecs( int size );
		void SolveFirstOrder( Field* velocity, Field* pressure, double nu, int svvCutoff );
		void SolveSecondOrder( Field* velocity, Field* pressure, Field* prevVel, double nu, int svvCutoff );
};
