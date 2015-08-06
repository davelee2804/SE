class PotentialVorticityEqn {
	public:
		PotentialVorticityEqn( Field* _omega, Field* _psi, Field* _velocity, Field* _prevOmega, Field* _prevVel );
		~PotentialVorticityEqn();
		void Solve( double dt );
		Field* omega;
		Field* psi;
		Field* velocity;
		Field* prevOmega;
		Field* prevVel;
		Advector* advector;
		Vec x;
		Vec b;
		Mat A;
		Vector* x0;
		Vector* x1;
		Vector* b0;
		Vector* b1;
		Matrix* A00;
		Matrix* A10;
		Matrix* A11;
};
