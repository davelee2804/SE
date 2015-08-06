class QuasiGeostrophicEqn {
	public:
		QuasiGeostrophicEqn( Field* _phi, Field* _omega, double _F, double _beta, double _nu, double _r, double _dt, Field* _topo  );
		~QuasiGeostrophicEqn();
		Field* 		omega;
		Field* 		phi;
		double		F;
		double		beta;
		double		nu;
		double		r;
		double 		dt;
		Field*		topo;
		void 		Solve( Field* phiPrev, Field* omegaPrev );
		void		FirstOrder();
		void 		SecondOrder( Field* phiPrev, Field* omegaPrev );
		Field*		GenVel( Field* psi );
	private:
		Field*		helmRhs1;
		Field*		helmRhs2;
		Field*		alpha;
		Field*		alphaHat;
		Field*		phiHat;
		Field*		omegaHat;
		HelmholtzEqn*	he1;
		HelmholtzEqn*	he2;
		void		GenPV( Field* pv, Field* psi, Field* vort );
		/* dissipdative terms stuff */
		bool		secondOrderAssembled;
		Vec		b;
		Vec		f;
		Mat		M;
		Vector*		qRHS;
		void		AssembleQMat( double a );
		void		SolveQ();
};

class DisRHS : public RHSOp {
	public:
		DisRHS( std::string _name, Mesh* mesh, double _constant, Field* _field, Field* _omega, double _r, double _nu );
		virtual ~DisRHS();
		virtual void AssembleElement( int el_i, double* rhs );
	private:
		Field*		omega;
		double		r;
		double		nu;
		double**	gW;
};
