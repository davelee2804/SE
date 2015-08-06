double CalcTimeStep( Field* velocity, double courantNo );

double CalcViscosity( Mesh* mesh );

void WriteXDMFHeader( int timeStep );
void WriteXDMFFooter( int timeStep );
void WriteXDMF( Field** fields, int nFields, int timeStep, double time, double dt );
void WriteXDMFTemporal( int nTimeSteps, int dumpEvery );

double FieldError( Field* numeric, Field* analytic, int dof, bool normalise );

void CreateVector( Field* field, Vec* v );
void CreateMatrix( Field* rowField, Field* colField, Mat* A );

//void CreateNodalToModalMatrix( Field* field, Mat* N2M );
