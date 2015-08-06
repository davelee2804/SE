class SchurComplement {
	public:
		SchurComplement( Mat _Kinv, Mat _G, Mat _D, Mat _C );
		~SchurComplement();
		void Solve( std::string prefix, Vec v, Vec p, Vec f, Vec h );
	private:
		Mat Kinv;
		Mat G;
		Mat D;
		Mat C;
		Mat DKinv;
		Mat DKinvG;
};
