class SVV : public Operator {
        public:
                SVV( std::string _name, Field* _rowField, Field* _colField, double _constant, int _mN );
                virtual ~SVV();
                int mN;
                virtual void AssembleElement( int el_i, double* A );
        private:
		Jacobi** mBasis;
		double* Lxx;
		double* Lyy;
		double* B;
		double* Q;
		double* W;
		double* Bt;
		double*	BtW;
		double* BtWB;
		double* BtWBinv;
		double* M;
		double* Minv;
		double* Dx;
		double* Dy;
		double* Dxt;
		double* Dyt;
		double* MinvQ;
		double* MinvQM;
		double* MinvQMDx;
		double* MinvQMDy;
		double* MinvQMDxt;
		double* MinvQMDyt;
		double* MinvQMDxtW;
		double* MinvQMDytW;
};
