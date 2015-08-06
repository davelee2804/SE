class NeumannBC : public RHSOp {
	public:
		NeumannBC( std::string _name, Mesh* _mesh, Field* _velocity, Field* _prevVel, Field* _velHat, double _nu );
		~NeumannBC();
		Field* velocity;
		Field* prevVel;
		Field* velHat;
		double nu;
		void AssembleElement( int el_i, double* rhs );
	private:
		double* weight;
		double* abcissa;
		double* Ni;
		double** dwdx;
		Field* vorticity;
		double* ShapeFuncs( double x );
};
