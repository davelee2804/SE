class NavierStokesSpectral {
	public:
		NavierStokesSpectral( Field* _velocity, Field* _pressure, double _nu );
		~NavierStokesSpectral();
		Field* 		velocity;
		Field* 		pressure;
		void Solve( double dt, Field* velPrev );
		void MapToArray( Field* field, double* array, int dof );
		void MapFromArray( Field* field, double* array, int dof );
	private:
		double 		nu;
		int		nx;
		int		ny;
		int		nTot;
		fftw_plan	uHat_forward;
		fftw_plan	vHat_forward;
		fftw_plan	u_backward;
		fftw_plan	v_backward;
		fftw_plan	p_backward;
		double*		p_real;
		double*		u_real;
		double*		v_real;
		double*		uHat_real;
		double*		vHat_real;
		fftw_complex*	p_fourier;
		fftw_complex*	u_fourier;
		fftw_complex*	v_fourier;
		fftw_complex*	uHat_fourier;
		fftw_complex*	vHat_fourier;
		double SVV( int kx, int ky, int cutoff );
		void SolveFirstOrder( double dt );
		void SolveSecondOrder( double dt, Field* velPrev );
};
