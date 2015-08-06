class RHSOp {
	public:
		RHSOp( std::string _name, Mesh* _mesh, double _constant, Field* _field );
		virtual ~RHSOp();
		std::string name;
		Mesh* mesh;
		double constant;
		Field* field;
		virtual void AssembleElement( int el_i, double* rhs ) = 0;
};

class FieldRHS : public RHSOp {
	public:
		FieldRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field );
		virtual ~FieldRHS();
		virtual void AssembleElement( int el_i, double* rhs );
};

class GradFieldRHS : public RHSOp {
	public:
		GradFieldRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field );
		virtual ~GradFieldRHS();
		virtual void AssembleElement( int el_i, double* rhs );
};

class DivVelRHS : public RHSOp {
	public:
		DivVelRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field );
		virtual ~DivVelRHS();
		virtual void AssembleElement( int el_i, double* rhs );
};

class GradDotFieldRHS : public RHSOp {
	public:
		GradDotFieldRHS( std::string _name, Mesh* _mesh, double _constnat, Field* _field );
		~GradDotFieldRHS();
		virtual void AssembleElement( int el_i, double* rhs );
};

class ConvectionRHS : public RHSOp {
        public:
                ConvectionRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field );
                virtual ~ConvectionRHS();
                virtual void AssembleElement( int el_i, double* rhs );
	private:
		double **du;
};

class StressTensorRHS : public RHSOp {
	public:
		StressTensorRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field );
		virtual ~StressTensorRHS();
		virtual void AssembleElement( int el_i, double* rhs );
	private:
		double **du;
};

class DivHeightVelRHS : public RHSOp {
	public:
		DivHeightVelRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field, Field* _velocity );
		virtual ~DivHeightVelRHS();
		Field* velocity;
		virtual void AssembleElement( int el_i, double* rhs );
};

class LaplacianRHS : public RHSOp {
	public:
		LaplacianRHS( std::string _name, Mesh* _mesh, double _constant, Field* _field );
		virtual ~LaplacianRHS();
		virtual void AssembleElement( int el_i, double* rhs );
};

class BiharmonicSurfaceRHS : public RHSOp {
	public:
		BiharmonicSurfaceRHS( std::string _name, Mesh* _mehs, double _constant, Field* _field );
		virtual ~BiharmonicSurfaceRHS();	
		double** g2;
		virtual void AssembleElement( int el_i, double* rhs );
};
