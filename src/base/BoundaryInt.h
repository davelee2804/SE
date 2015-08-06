class BoundaryInt : public RHSOp {
	public:
		BoundaryInt( std::string _name, Mesh* _mesh, Field* _field, Field* _fPrev );
		~BoundaryInt();
		Field* fPrev;
		void AssembleElement( int el_i, double* rhs );
	private:
		double* weight;
		double* abcissa;
};
