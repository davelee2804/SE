#include <cmath>
#include "Jacobi.h"

using namespace std;

Jacobi::Jacobi( int _n, double _a ) {
	int i;

	n = _n;
	a = _a;
	c = new double[n+1];
	for( i = 0; i <= n; i++ ) {
		c[i] = 0.0;
	}
}

Jacobi::~Jacobi() {
	delete[] c;
}

void Jacobi::Gen( Jacobi* m1, Jacobi* m2 ) {
	int i;
	double c1 = 2*n*(2*a + n)*(2*a + 2*n - 2);
	double c2 = (2*a + 2*n - 1)*(2*a + 2*n)*(2*a + 2*n - 2);
	double c3 = 2*(a + n - 1)*(a + n - 1)*(2*a + 2*n);

	c2 /= c1;
	c3 /= c1;

	for( i = 0; i <= m1->n; i++ ) {
		c[i+1] += c2*m1->c[i];
	}
	for( i = 0; i <= m2->n; i++ ) {
		c[i] -= c3*m2->c[i];
	}
}

double Jacobi::Eval( double x ) {
	int i;
	double p = 0.0;

	for( i = n; i >= 0; i -= 2 ) {
		p += c[i]*pow(x,i);
	}
	return p;
}

double Jacobi::EvalTotal( double x ) {
	int i;
	double p = 0.0;

	for( i = n; i >= 0; i-- ) {
		p += c[i]*pow(x,i);
	}
	return p;
}

double Jacobi::EvalDeriv( double x ) {
	int i;
	double dp = 0.0;

	for( i = n; i >= 1; i -= 2 ) {
		dp += c[i]*i*pow(x,i-1);
	}
	return dp;
}

double Jacobi::Eval2ndDeriv( double x ) {
	int i;
	double d2p = 0.0;

	for( i = n; i >= 2; i -= 2 ) {
		d2p += c[i]*i*(i-1)*pow(x,i-2);
	}
	return d2p;
}
