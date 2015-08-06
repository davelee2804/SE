class TwoLayerSW_Energetics {
	public:
		TwoLayerSW_Energetics( TwoLayerSW* _sw );
		~TwoLayerSW_Energetics();
		double 	CalcKE( int ts, int layer );
		double 	CalcPE( int ts );
		double 	CalcPower_SurfPres( int ts, int layer );
		double 	CalcPower_WindStress( int ts );
		double 	CalcPower_Friction( int ts );
		double 	CalcPower_Viscosity( int ts, int layer );
		void	Write( int ts );
	private:
		TwoLayerSW*	sw;
		void		CalcHeights( TwoLayerParams* params, Field* eta, Field** height1, Field** height2 );
};
