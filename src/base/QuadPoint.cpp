#include <cstdlib>
#include "QuadPoint.h"

QuadPoint::QuadPoint() {
	coord = NULL;
	weight = 0.0;
}

QuadPoint::QuadPoint( double* _coord, double _weight ) {
	coord = _coord;
	weight = _weight;
}

QuadPoint::~QuadPoint() {
	delete[] coord;
}
