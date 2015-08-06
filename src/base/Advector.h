class Advector {
	public:
		Advector( Field* _field, Field* _velocity );
		Advector( Field* _field, Field* _velocity, Field* _fieldMinusOne, Field* _velMinusOne );
		~Advector();
		void Advect( double dt );
		void Interpolate( Field* field, double* coord, double* psi );
		Field* field;
		Field* velocity;
		Field* fieldMinusOne;
		Field* velMinusOne;
		Field* velPlusOne;
		Field* velPlusHalf;
		Field* velMinusHalf;
		Field* fieldSL;
		bool secondOrder;
		bool cubicInterp;
	private:
		double** ptsX;
		double** ptsY;
		void PeriodicUpdate( double* coord );
		void IntegrateRK2( double dt, Field* velField, double* origin, double* coord );
		void IntegrateRK4( double dt, double* origin, double* coord );
		void FirstOrder( double dt );
		void SecondOrder( double dt );
		void InterpLagrange( double x, double* coords, int nDofs, double** vals, double* psi );
		void ExtrapolateVelocities();
};
