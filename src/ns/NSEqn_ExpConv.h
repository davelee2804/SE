class NSEqn_ExpConv {
	public:
		NSEqn_ExpConv( Field* velocity, Field* pressure, double _dt );
		~NSEqn_ExpConv();
		Field* velStep1;
		Field* velStep2;
		double dt;
		Mat    Av;
		Mat    Ap;
		Vec    x;
		Vec    b;
		Vec    bv;
		Vec    bp;
		bool   firstOrderInit;
		void Solve( Field* velocity, Field* pressure, Field* prevVel, bool doVelMat, double nu, int svvCutoff );
	private:
		void InitMats( int vSize, int pSize, int nVelNodes, int nPresNodes );
		void InitVecs( int size );
		void SolveFirstOrder( Field* velocity, Field* pressure, double nu, int svvCutoff );
		void SolveSecondOrder( Field* velocity, Field* pressure, Field* prevVel, bool doVelMat, double nu, int svvCutoff );
};
