class Energetics_SW {
	public:
		Energetics_SW( ShallowWaterEqn* _sw );
		~Energetics_SW();
		double 	CalcKE( int ts );
		double 	CalcPE( int ts );
		double 	CalcPower_WindStress( int ts );
		double 	CalcPower_Friction( int ts );
		double 	CalcPower_Viscosity( int ts );
		void	Write( int ts );
	private:
		ShallowWaterEqn*	sw;
};
