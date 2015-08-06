class QuadPoint {
	public:
		QuadPoint();
		QuadPoint( double* _coord, double _weight );
		~QuadPoint();
		double* coord;
		double weight;
};
