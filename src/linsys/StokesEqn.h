class StokesEqn {
	public:
		StokesEqn( Field* _velocity, Field* _pressure, Field* _forcing );
		~StokesEqn();
		Field* velocity;
		Field* pressure;
		Field* forcing;
		Vec x;
		Vec b;
		Mat A;
		Mat P;
		Vector* v;
		Vector* p;
		Vector* f;
		Vector* h;
		Matrix* K;
		Matrix* G;
		Matrix* D;
		Matrix* C;
		void Solve();
};
