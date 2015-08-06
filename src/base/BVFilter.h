class BVFilter {
	public:
		BVFilter( Field* _field, double _nu, int _p, int _s );
		~BVFilter();
		Field*  	nodalField;		/* the field to be filtered */
		Field*		modalField;		/* the field in modal space */
		double  	nu;			/* the filter 'viscosity' */
		int     	p;			/* controls the steepness of the filter */
		int		s;			/* modal frequency at which the filter kicks in */
		void		Apply();
	private:
		double*		B;
		double*		Binv;
		double*		w_k;
		double		Sigma( int k );
};
