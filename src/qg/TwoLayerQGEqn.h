class TwoLayerQGEqn {
	public:
		TwoLayerQGEqn( Field* _omega1, Field* _psi1, Field* _omega2, Field* _psi2, double _F1, double _F2, double _beta, double _E2, double _r, double _nu, Field* tau_e, double _dt );
		~TwoLayerQGEqn();
		Field*		omega1;
		Field*		omega2;
		Field*		psi1;
		Field*		psi2;
		double		F1;
		double		F2;
		double		beta;
		double		E2;
		double		r;
		double		nu;
		double		dt;
		void		Solve( Field* phi1Prev, Field* omega1Prev, Field* phi2Prev, Field* omega2Prev );
		Field*		GenVel( Field* psi, std::string name );
	private:
		Field*		tau1;		// baroclinic stream function: psi1 - psi2
		Field*		tau2;		// baroclinic stream function: psi1 - psi2
		Field*		chi1;		// barotropic stream function: psi1 + psi2
		Field*		chi2;		// barotropic stream function: psi1 + psi2
		Field*		omega1Hat;
		Field*		omega2Hat;
		Field*		psi1Hat;
		Field*		psi2Hat;
		Field*		helmRhs1;
		Field*		helmRhs2;
		Field*		fishRhs1;
		Field*		fishRhs2;
		HelmholtzEqn*	helmEqn1;
		HelmholtzEqn*	helmEqn2;
		PoissonEqn*	fishEqn1;
		PoissonEqn*	fishEqn2;
		Field*		alpha1;
		Field*		alpha2;
		Field*		alpha1Hat;
		Field*		alpha2Hat;
		Field*		phi1Hat;
		Field*		phi2Hat;
		Mat		M1;
		Mat		M2;
		Vec		b1;
		Vec		b2;
		Vec		f1;
		Vec 		f2;
		Vector*		q1RHS;
		Vector*		q2RHS;
		bool		secondOrderAssembled;
		void		GenPV( Field* q, Field* omega, Field* psi_i, Field* psi_j, double F );
		void		AssembleQMat( Mat* M, Vec* b, Field* qi, double a );
		void		SolveQ( Mat M, Vec f, Vec b, Vector* rhs, Field* qi );
		void		FirstOrder();
		void		SecondOrder( Field* phi1Prev, Field* omega1Prev, Field* phi2Prev, Field* omega2Prev );
};

class DissipationRHS : public RHSOp {
	public:
		DissipationRHS( std::string _name, Mesh* mesh, double _constant, Field* _field, Field* _tau, Field* _omegaHat, double _rF, double _Ei, Field* _tau_e, double _nu );
		virtual ~DissipationRHS();
		virtual void AssembleElement( int el_i, double* rhs );
	private:
		Field*		tau;
		Field* 		omegaHat;
		double		rF;
		double		Ei;
		Field*		tau_e;
		double		nu;
		double**	gW;
};
