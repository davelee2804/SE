class QuasiGeostrophicSpectral {
	public:
		QuasiGeostrophicSpectral( Field* _vorticity, Field* _streamFunc, Field* _velocity, double _beta );
		~QuasiGeostrophicSpectral();
		Field* 		vorticity;
		Field* 		streamFunc;
		Field*		velocity;
		double		beta;
		void Solve( double dt, Field* vortPrev, Field* velPrev );
		void MapToArray( Field* field, double* array, int dof );
		void MapFromArray( Field* field, double* array, int dof );
	private:
		double 		nu;
		int		nx;
		int		ny;
		int		nTot;
		fftw_plan	wHat_forward;
		fftw_plan	p_backward;
		fftw_plan	w_backward;
		fftw_plan	u_backward;
		fftw_plan	v_backward;
		double*		p_real;
		double*		w_real;
		double*		u_real;
		double*		v_real;
		double*		wHat_real;
		fftw_complex*	p_fourier;
		fftw_complex*	w_fourier;
		fftw_complex*	u_fourier;
		fftw_complex*	v_fourier;
		fftw_complex*	wHat_fourier;
};
