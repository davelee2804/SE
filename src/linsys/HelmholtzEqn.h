class HelmholtzEqn {
	public:
		HelmholtzEqn( Field* _solField, Field* _rhsField, double _c, double _nu, double _dt, int _svvCutoff );
		~HelmholtzEqn();
		Field* solField;
		Field* rhsField;
		double c;
		double nu;
		double dt;
		int svvCutoff;
		Vector* sol;
		Vector* rhs;
		Matrix* op;
		Vec x;
		Vec f;
		Vec b;
		Mat A;
		void Assemble();
		void Solve( std::string name, bool removeNullSpace );
	private:
		bool assembled;
};
