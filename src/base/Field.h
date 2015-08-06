typedef double Coord[3];

typedef double ( BCFunc ) ( Coord coord );

class Field {
	public:
		Field( std::string _name, Mesh* _mesh, int _nDofs, BCs* _bcs );
		~Field();
		std::string name;
		Mesh* mesh;
		double** vals;
		int nDofs;
		BCs* bcs;
		double** GNix;
		void   Copy( Field* field );
		void   InterpLocal( int el, Coord lCoord, double* val );
		void   InterpGlobal( Coord gCoord, double* val );
		void   InterpDerivsLocal( int el, Coord lCoord, double** gradVal );
		void   InterpDerivsGlobal( Coord gCoord, double** gradVal );
		void   InterpDerivsWithGlobalShapeFuncDerivs( int el, double** globalGNix, double** gradVal );
		void   InterpSecondDerivs( int el_i, int pt_i, double** g2Val );
		void   InterpSecondDerivsAtCoord( int el_i, double* coord, double** g2Val );
		double Integrate( int dof, bool average );
		void   SetBCConst( std::string side, int dof, double val );
		void   SetBCFunc( std::string side, int dof, BCFunc* bcFunc );
		void   SetICConst( int dof, double val );
		void   SetICFunc( int dof, BCFunc* icFunc );
		void   PeriodicUpdate();
		void   Save( int timeStep );
		void   Read( int timeStep );
};

