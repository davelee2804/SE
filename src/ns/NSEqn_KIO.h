class NSEqn_KIO {
	public:
		NSEqn_KIO( Field* _velocity, Field* _pressure, double _dt, double _nu );
		~NSEqn_KIO();
		Field* 		velocity;
		Field* 		pressure;
		double 		dt;
		double		nu;
		void 		Setup( int order );
		void 		Solve( Field* velPrev );
	private:
		Mat	M_u;
		Mat	M_p;
		Vec	b_u;
};
