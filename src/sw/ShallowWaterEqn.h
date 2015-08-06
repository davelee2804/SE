/* Note: All parameters are assumed to have been divided by Ro (the Rossby number) = U/(f_0*L)
 *       as the advective term is scaled to unity */

typedef struct {
	double L;	/* horizontal length scale	(m) 		*/
        double H;	/* vertical length scale	(m) 		*/
        double U;	/* velocity scale 		(m/s)		*/
        double f0;	/* coriolis scale 		(1/s)		*/
        double beta;	/* coriolis gradient 		(1/ms)		*/
        double g;	/* gravity			(1/ms^2) 	*/
        double gamma;	/* friction coefficient		(1/s) 		*/
        double tau;	/* wind stress coefficient	(kg/ms^2) 	*/
	double kws;	/* wind stress vertical period	(1/m)		*/
        double nu;	/* viscosity coefficient 	(m^2/s)		*/
        double rho;	/* density			(kg/m^3) 	*/
} SWParams;

class ShallowWaterEqn {
	public:
		ShallowWaterEqn( Field* _velocity, Field* _pressure, double _dt, int _svvCutoff, Field* _topography, SWParams* _params );
		~ShallowWaterEqn();
		Field* 		velocity;
		Field* 		pressure;
		SWParams* 	params;
		double 		dt;
		int    		svvCutoff;
		double		qgScale;
		Field* 		topography;
		bool 		uNullSp;
		bool 		vNullSp;
		bool 		pNullSp;
		bool		presDirichlet;
		void 		Solve( Field* velPrev, Field* presPrev );
		void 		SecondOrder( Field* velPrev, Field* presPrev );
		void		FirstOrder();
		void 		Assemble( int order );
		Mat 		A;
		Vec 		b;
		Vec 		x;
		Vec 		f;
	private:
		void 		InitObjs();
		void		SolveLinAlg( Vec f, Vec x );
};
