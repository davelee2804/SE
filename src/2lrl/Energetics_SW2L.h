class Energetics_SW2L {
	public:
		Energetics_SW2L( SW2L* _sw, BarotropicEqn* _bt );
		~Energetics_SW2L();
		double 	CalcKE( int ts, int layer );
		double 	CalcPE( int ts );
		double 	CalcPower_SurfPres( int ts, int layer );
		double 	CalcPower_WindStress( int ts );
		double 	CalcPower_Friction( int ts );
		double 	CalcPower_Viscosity( int ts, int layer );
		void	Write( int ts );
	private:
		SW2L*		sw;
		BarotropicEqn*	bt;
};
