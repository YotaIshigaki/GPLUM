#ifndef BSPLINE_H
#define BSPLINE_H
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <cassert>

namespace PMEModule {

struct BSpline {
    public:
	typedef BSpline* Ptr;

	explicit BSpline(int _n): n(_n), f(n), df(n) {
	  if (n < 2) {
	    throw("order is too small");
	  }
	}
	~BSpline() {}

	void calculate_nodiff(double u)
	{
	    assert(n >= 2);
	    double x = u - floor(u);
	    f[0] = x;
	    f[1] = 1.0-x;
	    for (int i=2;i<n;++i) {
		recurrence(x, i);
	    }
	}
	void calculate(double u)
	{
	    assert(n >= 2);
	    double x = u - floor(u);
	    f[0] = x;
	    f[1] = 1.0-x;
	    df[0] = 1.0;
	    df[1] = -1.0;
	    for (int i=2;i<n-1;++i) {
		recurrence(x, i);
	    }
	    difference();
	    recurrence(x, n-1);
	}
	double function(double u)
	{
	    calculate(u);
	    if (u < 0 || u >= n) {
		return 0.0;
	    }
	    return f[int(floor(u))];
	}
	double derivative(double u)
	{
	    calculate(u);
	    if (u < 0 || u >= n) {
		return 0.0;
	    }
	    return df[int(floor(u))];
	}
	double* getFunction() { return &f[0]; }
	double* getDerivative() { return &df[0]; }
    private:
	int n;
	std::vector<double> f;
	std::vector<double> df;
	void recurrence(double x, int i)
	{
	    f[i] = (1-x)*f[i-1]/i;
	    for (int j=i-1;j>0;--j) {
		f[j] = ((x+j)*f[j]+(i-j+1-x)*f[j-1])/i; //  Ulrich Essmann 1995 (4.1)
	    }
	    f[0] = x*f[0]/i;
	}
	void difference() //  Ulrich Essmann 1995 (4.2)
	{
	    df[0] = f[0];
	    for (int i=1;i<n-1;++i) {
		df[i] = f[i]-f[i-1];
	    }
	    df[n-1] = -f[n-2];
	}
};
}
#endif
