#pragma once
/* C++ headers */
#include "common.h"
/* User-defined headers */

// [Notes]
// (1) The best answer at this link (https://oshiete.goo.ne.jp/qa/8674219.html)
//     gives a good description about the usage of std::function<>.


template <class real_t>
void polint(real_t xa[], real_t ya[],
            const int n, const real_t x,
            real_t & y, real_t & dy) {
    int ns = 0;
    real_t dif = std::fabs(x-xa[0]);
    std::vector<real_t> c(n), d(n);
    for (int i = 0; i < n; i++) {
        const real_t dift = std::fabs(x-xa[i]);
        if (dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns--];
    for (int m = 1; m < n; m++) {
        for (int i = 0 ; i < n-m; i++) {
            const real_t ho = xa[i] - x;
            const real_t hp = xa[i+m] - x;
            const real_t w = c[i+1] - d[i];
            real_t den = ho - hp;
            if (den == 0.0) {
                std::cerr << "Error in routine polint" << std::endl;
                assert(false);
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
    }
}

template <class real_t>
real_t trapzd(std::function<real_t(const real_t)> func,
              const real_t a, const real_t b, const int n) {
    static real_t s;

    if (n == 1) {
        return (s=0.5*(b-a)*(func(a)+func(b)));
    } else {
        int it {1};
        for (int j = 1; j < n-1; j++) it <<= 1;
        const real_t tnm = it;
        const real_t del = (b-a)/tnm;
        real_t x = a + 0.5*del;
        real_t sum {0.0};
        for (int j = 0; j < it; j++, x += del) sum += func(x);
        s = 0.5*(s+(b-a)*sum/tnm); 
        return s;
    }
}

template <class real_t>
real_t qromb(std::function<real_t(const real_t)> func, const real_t a, const real_t b) {
    constexpr real_t EPS = 1.0e-6;
    constexpr int JMAX = 20;
    constexpr int JMAXP = (JMAX+1); 
    constexpr int K = 5;
    real_t ss,dss;
    real_t s[JMAX], h[JMAXP];
    h[0]=1.0;
    for (int j = 1; j <= JMAX; j++) {
        s[j-1] = trapzd(func, a, b, j);
        if (j >= K) {
            polint(&h[j-K], &s[j-K], K, 0.0, ss, dss);
            if (std::fabs(dss) <= EPS*std::fabs(ss)) return ss;
        }
        h[j]=0.25*h[j-1];
    }
    std::cerr << "Too many steps in routine qromb" << std::endl;
    return 0.0;
}

template <class real_t>
void gauleg(const real_t x1, const real_t x2, real_t x[], real_t w[], const int n) {
    constexpr real_t EPS = 1.0e-14;
    constexpr real_t PI = 3.141592653589793238462643383279;
    const int m=(n+1)/2;
    const real_t xm=0.5*(x2+x1);
    const real_t xl=0.5*(x2-x1);
    for (int i=0; i<m; i++) {
        real_t z1,pp;
        real_t z = std::cos(PI*(i+0.75)/(n+0.5));
        do {
            real_t p1 = 1.0;
            real_t p2 = 0.0;
            for (int j=0; j<n; j++) {
                real_t p3 = p2;
                p2 = p1;
                p1 = ((2.0*j+1.0)*z*p2-j*p3)/(j+1);
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (std::fabs(z-z1) > EPS);
        x[i] = xm-xl*z;
        x[n-1-i] = xm+xl*z;
        w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
        w[n-1-i] = w[i];
    }
}

template <class real_t>
real_t qgaus(std::function<real_t(const real_t)> func, const real_t a, const real_t b) {
    static real_t x[] = {0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
    static real_t w[] = {0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};
    const real_t xm = 0.5*(b+a);
    const real_t xr = 0.5*(b-a);
    real_t s {0.0}; 
    for (int j = 1; j <= 5; j++) {
        const real_t dx = xr * x[j];
        s += w[j]*(func(xm+dx) + func(xm-dx));
    }
    return s *= xr;
}

template <class real_t, int jmax>
real_t qgausap(std::function<real_t(const real_t)> func, const real_t a, const real_t b) {
    static_assert(jmax > 0, "jmax must be a positive integer.");
    static bool is_initialized {false};
    static std::vector<real_t> x(jmax+1), w(jmax+1);
    if (is_initialized == false) {
        std::vector<real_t> x_tmp(2*jmax), w_tmp(2*jmax);
        gauleg<real_t>(-1.0, 1.0, &x_tmp[0], &w_tmp[0], 2*jmax);
        x[0] = 0.0; w[0] = 0.0;
        for (int i = 1; i <= jmax; i++) {
            x[i] = x_tmp[i + jmax - 1];
            w[i] = w_tmp[i + jmax - 1];
        }
    }
    const real_t xm = 0.5*(b+a);
    const real_t xr = 0.5*(b-a);
    real_t s {0.0}; 
    for (int j = 1; j <= 5; j++) {
        const real_t dx = xr * x[j];
        s += w[j]*(func(xm+dx) + func(xm-dx));
    }
    return s *= xr;
}

template <class real_t>
class quad2d {
private:
    real_t xsav_;
    std::function<real_t(const real_t, const real_t)> func_;
    std::function<real_t(const real_t)> y1_;
    std::function<real_t(const real_t)> y2_;
    real_t (*intg_)(std::function<real_t(const real_t)>, const real_t, const real_t);
    
    real_t f1(const real_t x) {
        this->xsav_ = x;
        return this->intg_([this](const real_t y) -> real_t { return this->f2(y);},
                           this->y1_(x),
                           this->y2_(x));
    }

    real_t f2(const real_t y) {
        return func_(this->xsav_, y); 
    }

public:

    real_t operator () (std::function<real_t(const real_t, const real_t)> func,
                        const real_t x1,
                        const real_t x2,
                        std::function<real_t(const real_t)> y1,
                        std::function<real_t(const real_t)> y2,
                        real_t (*intg)(std::function<real_t(const real_t)>, const real_t, const real_t)) {
        this->func_ = func;
        this->y1_ = y1;
        this->y2_ = y2;
        this->intg_ = intg;
        return this->intg_([this](const real_t x) -> real_t { return this->f1(x);},
                           x1, x2);
    }

    void test() {
        const real_t rc = 5.0;
        const real_t val_qgaus
            = (*this)([](const real_t x, const real_t y)
                        -> real_t {return std::exp(-(x*x+y*y));},
                      -rc, rc,
                      [](const real_t x) -> real_t {return -std::sqrt(rc*rc - x*x);},
                      [](const real_t x) -> real_t {return  std::sqrt(rc*rc - x*x);},
                      qgaus<real_t>);
        const real_t val_qgausap
            = (*this)([](const real_t x, const real_t y) 
                        -> real_t {return std::exp(-(x*x+y*y));},
                      -rc, rc,
                      [](const real_t x) -> real_t {return -std::sqrt(rc*rc - x*x);},
                      [](const real_t x) -> real_t {return  std::sqrt(rc*rc - x*x);},
                      qgausap<real_t,30>);
        std::cout << "Calculated results: " << std::endl;
        std::cout << "     qgaus = " << val_qgaus << std::endl;
        std::cout << "   qgausap = " << val_qgausap << std::endl;
        std::cout << "Exact sol. = " << 3.1416 << std::endl;
    }

};


template <class real_t>
class quad3d {
private:
    real_t xsav_, ysav_;
    std::function<real_t(const real_t, const real_t, const real_t)> func_;
    std::function<real_t(const real_t)> y1_;
    std::function<real_t(const real_t)> y2_;
    std::function<real_t(const real_t, const real_t)> z1_;
    std::function<real_t(const real_t, const real_t)> z2_;
    real_t (*intg_)(std::function<real_t(const real_t)>, const real_t, const real_t);
    
    real_t f1(const real_t x) {
        this->xsav_ = x;
        return this->intg_([this](const real_t y) -> real_t { return this->f2(y);},
                           this->y1_(x),
                           this->y2_(x));
    }

    real_t f2(const real_t y) {
        this->ysav_ = y;
        return this->intg_([this](const real_t z) -> real_t { return this->f3(z);},
                           this->z1_(this->xsav_, y),
                           this->z2_(this->xsav_, y)); 
    }

    real_t f3(const real_t z) {
        return func_(this->xsav_, this->ysav_, z); 
    }

public:

    real_t operator () (std::function<real_t(const real_t, const real_t, const real_t)> func,
                        const real_t x1,
                        const real_t x2,
                        std::function<real_t(const real_t)> y1,
                        std::function<real_t(const real_t)> y2,
                        std::function<real_t(const real_t, const real_t)> z1,
                        std::function<real_t(const real_t, const real_t)> z2,
                        real_t (*intg)(std::function<real_t(const real_t)>, const real_t, const real_t)) {
        this->func_ = func;
        this->y1_ = y1;
        this->y2_ = y2;
        this->z1_ = z1;
        this->z2_ = z2;
        this->intg_ = intg;
        return this->intg_([this](const real_t x) -> real_t { return this->f1(x);},
                           x1, x2);
    }

    void test() {
        const real_t rc = 5.0;
        const real_t val_qgaus
            = (*this)([](const real_t x, const real_t y, const real_t z) -> real_t {return std::exp(-(x*x+y*y+z*z));},
                      -rc, rc,
                      [](const real_t x) -> real_t {return -std::sqrt(rc*rc - x*x);},
                      [](const real_t x) -> real_t {return  std::sqrt(rc*rc - x*x);},
                      [](const real_t x, const real_t y) -> real_t {return -std::sqrt(rc*rc - x*x - y*y);},
                      [](const real_t x, const real_t y) -> real_t {return  std::sqrt(rc*rc - x*x - y*y);},
                      qgaus<real_t>);
        const real_t val_qgausap
            = (*this)([](const real_t x, const real_t y, const real_t z) -> real_t {return std::exp(-(x*x+y*y+z*z));},
                      -rc, rc,
                      [](const real_t x) -> real_t {return -std::sqrt(rc*rc - x*x);},
                      [](const real_t x) -> real_t {return  std::sqrt(rc*rc - x*x);},
                      [](const real_t x, const real_t y) -> real_t {return -std::sqrt(rc*rc - x*x - y*y);},
                      [](const real_t x, const real_t y) -> real_t {return  std::sqrt(rc*rc - x*x - y*y);},
                      qgausap<real_t,15>);
        std::cout << "Calculated results: " << std::endl;
        std::cout << "     qgaus = " << val_qgaus << std::endl;
        std::cout << "   qgausap = " << val_qgausap << std::endl;
        std::cout << "Exact sol. = " << 5.56833 << std::endl;
    }

};
