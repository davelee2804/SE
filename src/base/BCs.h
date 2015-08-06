class BCs {
	public:
		BCs( bool _bcBottom, bool _bcTop, bool _bcLeft, bool _bcRight, Mesh* _mesh, int _vecOffset );
		BCs( bool* _bcBottom, bool* _bcTop, bool* _bcLeft, bool* _bcRight, Mesh* _mesh, int _vecOffset );
		~BCs();
		Mesh* mesh;
		int vecOffset;
		int* vecToFieldMap;
		int* fieldToVecMap;
		int size[2];
		bool bcBottom[2];
		bool bcTop[2];
		bool bcLeft[2];
		bool bcRight[2];
		int* bottom;
		int* top;
		int* left;
		int* right;
		bool IsBCNode( int node );
		bool IsBCNodeAndDof( int node, int dof );
		int* GetSide( std::string side, int* nNodes );
};
