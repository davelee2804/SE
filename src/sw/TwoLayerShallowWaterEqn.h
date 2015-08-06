/* Note: All parameters are assumed to have been divided by Ro (the Rossby number) = U/(f_0*L)
 *       as the advective term is scaled to unity */

class TwoLayerShallowWaterEqn {
	public:
		TwoLayerShallowWaterEqn( Field* _velocity1, 
					 Field* _pressure1, 
					 Field* _velocity2,
					 Field* _pressure2,
					 double _dt,
					 int    _svvCutoff,
					 Field* _topography,
					 SWParams* params,
					 double _gPrime );
		~TwoLayerShallowWaterEqn();
		Field* velocity1;
		Field* pressure1;
		Field* velocity2;
		Field* pressure2;
		double L;
		double H;
		double U;
		double f0;
		double beta;
		double g;
		double gPrime;
		double gamma;
		double tau;
		double kws;
		double nu;
		double rho;
		double dt;
		int    svvCutoff;
		Field* topography;
		bool uNullSp;
		bool vNullSp;
		bool pNullSp;
		int nNonLinIts;
		void Solve( Field* velPrev1, Field* presPrev1, Field* velPrev2, Field* presPrev2, int assMatOrder );
		void Assemble( double a );
		void CalcEnergetics( int timeStep, ShallowWaterEnergetics* topEnergetics, ShallowWaterEnergetics* bottomEnergetics );
	private:
		Mat A;
		Vec b;
		Vec x;
		Vec f;
		int presVecOffset1;
		int presVecOffset2;
		void InitObjs();
		void _Solve( Field* velPrev1, Field* presPrev1, Field* velPrev2, Field* presPrev2 );
};
