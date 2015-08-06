class FourierChebyshev {
	public:
		FourierChebyshev( Field* _inField, int _dof, int* _nx );
		~FourierChebyshev();
		void 		Trans();
		void		Write( std::string filename, double minAmp );
		Mesh*		mesh;
		Field*		field;
	private:
		ChebTrans*	ct;
		Field*		inField;
		int		dof;
		int		nx[2];
		double*		grid;
		fftw_complex*	tran;
		double*		real_y;
		double*		real_x;
		fftw_complex*	four_x;
		fftw_plan	fft_x;
		void		ModifyMesh();
		void		InterpField();
		void		FieldToGrid();
		void		GridTranspose( double* g, double* gt, int n, int m );
		void		GetRow( double* g, double* r, int row, int n );
		void		SetRow( double* g, double* r, int row, int n );
		void		SetRowComplex( fftw_complex* g, fftw_complex* r, int row, int n );
		void		Norm();
};
