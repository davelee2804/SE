/* Note: All parameters are assumed to have been divided by Ro (the Rossby number) = U/(f_0*L)
 *       as the advective term is scaled to unity */

class LinearShallowWaterEqn {
	public:
		LinearShallowWaterEqn( Field* _velocity, 
				 Field* _pressure, 
				 double _f0, 
				 double _beta, 
				 double _g, 
				 double _gamma, 
				 double _tau, 
				 double _nu,
				 double _dt, 
                                 int    _svvCutoff,
				 Field* _topography );
		~LinearShallowWaterEqn();
		Field* velocity;
		Field* pressure;
		double f0;
		double beta;
		double g;
		double gamma;
		double tau;
		double nu;
		double dt;
		int    svvCutoff;
		Field* topography;
		bool uNullSp;
		bool vNullSp;
		bool pNullSp;
		void Solve( Field* velPrev, Field* presPrev, int assMatOrder );
		void CalcEnergetics( int timeStep, double* KE, double* PE, double* visc, double* fric, double* wind );
	private:
		Mat A;
		Vec b;
		Vec x;
		Vec f;
		void _Solve( Field* velPrev, Field* presPrev );
		void InitObjs();
		void Assemble( double a );
};

