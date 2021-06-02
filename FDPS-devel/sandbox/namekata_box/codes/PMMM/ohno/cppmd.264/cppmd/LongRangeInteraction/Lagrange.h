#ifndef LAGRANGE_H
#define LAGRANGE_H
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>

namespace PMEModule {

struct Lagrange {
    public:
	typedef Lagrange* Ptr;

	Lagrange(unsigned int _p): p(_p), f(2*p), df(2*p) {}
	~Lagrange() {}

	void calculate_nodiff(double u)
	{
	    double x = u-floor(u)-p;
	    int n = 2*p;
	    double *q = &f[0];
	    for (int k=0;k<n;++k,++q,++x) {
		*q = 1.0;
		for (int j=1;j<n-k;++j) {
		    *q *= (x+j)/j;
		}
		for (int j=1;j<=k;++j) {
		    *q *= (j-x)/j;
		}
	    }
	}
	void calculate(double u)
	{
	    double x = u-floor(u)-p;
	    int n = 2*p;
	    double *q = &f[0];
	    double *dq = &df[0];
	    for (int k=0;k<n;++k,++q,++dq,++x) {
		*q = 1.0;
		*dq = 0.0;
		for (int j=1;j<n-k;++j) {
		    *dq = ((x+j)*(*dq)+(*q))/j;
		    *q *= (x+j)/j;
		}
		for (int j=1;j<=k;++j) {
		    *dq = ((j-x)*(*dq)-(*q))/j;
		    *q *= (j-x)/j;
		}
	    }
	}
	double function(double u)
	{
	    if (u < -p || u >= p) {
		return 0.0;
	    }
	    calculate(u);
	    return f[int(floor(u))+p];
	}
	double derivative(double u)
	{
	    calculate(u);
	    if (u < -p || u >= p) {
		return 0.0;
	    }
	    return df[int(floor(u))+p];
	}
	double* getFunction() { return &f[0]; }
	double* getDerivative() { return &df[0]; }
    private:
	int p;
	std::vector<double> f;
	std::vector<double> df;
};
}
#endif
