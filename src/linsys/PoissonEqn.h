class PoissonEqn {
	public:
		PoissonEqn( Field* _solField, Field* _rhsField, double _c );
		~PoissonEqn();
		Field* solField;
		Field* rhsField;
		double c;
		Vector* x0;
		Vector* f0;
		Matrix* A00;
		Vec x;
		Vec f;
		Vec b;
		Mat A;
		void Assemble();
		void Solve( bool removeNullSpace );
	private:
		bool assembled;
};
