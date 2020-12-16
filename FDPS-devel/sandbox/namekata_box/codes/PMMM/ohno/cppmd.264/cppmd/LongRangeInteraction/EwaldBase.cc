#include <iostream>
#include "EwaldBase.h"
#ifdef ICL
#include "MKLErfc.h"
#endif

using namespace std;
using namespace EwaldModule;

namespace { // This part is imported from ks_math.c. It must be rewritten by using C++.

const double KS_PI_PI = 9.86960440108935861883;
const double NR_GOLD = 1.618043;
const double NR_GLIMIT = 100.0;
const double NR_TINY = 1.0e-20;
const double NR_GOLD_R = 0.61803399;
const double NR_GOLD_C = (1.0-NR_GOLD_R);

#define NR_SHIFT2(a,b,c) (a)=(b);(b)=(c);
#define NR_SHIFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
static double maxarg1,maxarg2;
#define NR_MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1)>(maxarg2)?(maxarg1):(maxarg2))
#define NR_SIGN(a,b) ((b) >= 0.0 ? fabs(a):-fabs(a))

void ks_nr_mnbrake(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
		    double (*func)(double))
{
  double ulim,u,r,q,fu,dum;

  *fa = (*func)(*ax);
  *fb = (*func)(*bx);

  if(*fb > *fa){
    NR_SHIFT3(dum,*ax,*bx,dum);
    NR_SHIFT3(dum,*fa,*fb,dum);
  }
  *cx = (*bx)+NR_GOLD*(*bx-*ax);
  *fc = (*func)(*cx);
  /*
  printf("a %f b %f c %f  fa %f fb %f fc %f %f:%f\n",*ax,*bx,*cx,*fa,*fb,*fc
	 ,(*bx-*ax)/(*cx-*ax)
	 ,(*cx-*bx)/(*cx-*ax)
	 );
  */
  while(*fb > *fc){
    r = (*bx-*ax)*(*fb-*fc);
    q = (*bx-*cx)*(*fb-*fa);
    u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*NR_SIGN(NR_MAX(fabs(q-r),NR_TINY),q-r));
    /*    printf("u %f\n",u);*/
    ulim = (*bx)+NR_GLIMIT*(*cx-*bx);
    if((*bx-u)*(u-*cx) > 0.0){
      /*      printf("case 1\n");*/
      fu = (*func)(u);
      if(fu < *fc){
	*ax = (*bx);
	*bx = u;
	*fa = (*fb);
	*fb = fu;
	return;
      } else if(fu > *fb){
	*cx = u;
	*fc = fu;
	return;
      }
      u = (*cx)+NR_GOLD*(*cx-*bx);
      fu = (*func)(u);
    } else if((*cx-u)*(u-ulim) > 0.0){
      /*      printf("case 2\n");*/
      fu = (*func)(u);
      if(fu < *fc){
	NR_SHIFT3(*bx,*cx,u,*cx+NR_GOLD*(*cx-*bx));
	NR_SHIFT3(*fb,*fc,fu,(*func)(u));
      }
    } else if((u-ulim)*(ulim-*cx) >= 0.0){
      /*      printf("case 3\n");*/
      u = ulim;
      fu=(*func)(u);
    } else {
      /*      printf("case 4\n");*/
      u = (*cx)+NR_GOLD*(*cx-*bx);
      fu = (*func)(u);
    }
    NR_SHIFT3(*ax,*bx,*cx,u);
    NR_SHIFT3(*fa,*fb,*fc,fu);
  }
}
double ks_nr_golden(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin)
{
  double f1,f2,x0,x1,x2,x3;

  x0 = ax;
  x3 = cx;

  if(fabs(cx-bx) > fabs(bx-ax)){
    x1 = bx;
    x2 = bx+NR_GOLD_C*(cx-bx);
  } else {
    x2 = bx;
    x1 = bx-NR_GOLD_C*(bx-ax);
  }
  f1 = (*f)(x1);
  f2 = (*f)(x2);
  while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))){
    if(f2 < f1){
      NR_SHIFT3(x0,x1,x2,NR_GOLD_R*x1+NR_GOLD_C*x3);
      NR_SHIFT2(f1,f2,(*f)(x2));
    } else {
      NR_SHIFT3(x3,x2,x1,NR_GOLD_R*x2+NR_GOLD_C*x0);
      NR_SHIFT2(f2,f1,(*f)(x1));
    }
    /*    printf("%f %f  %f %f\n",x1,x2,f1,f2);*/
  }
  if(f1 < f2){
    *xmin = x1;
    return f1;
  } else {
    *xmin = x2;
    return f2;
  }
}
static double _solve_erfc_alpha;
static double _solve_erfc_cut_off;
static double _solve_erfc_sub;
static double _solve_erfc(double x)
{
  return fabs(erfc(_solve_erfc_alpha*x)/x-_solve_erfc_sub);
}
double ks_solve_erfc_cut_off(double alpha, double r)
{
  double a,b,c,fa,fb,fc;
  double min_x,min_y;

  a = 0.0001; 
  b = 1;
  _solve_erfc_alpha = alpha;
  _solve_erfc_sub = r;
  ks_nr_mnbrake(&a,&b,&c,&fa,&fb,&fc,_solve_erfc);
  min_y = ks_nr_golden(a,b,c,_solve_erfc,1e-5,&min_x);
  return min_x;
}
double _solve_erfc_alpha_func(double alpha)
{
  double a,b,c,fa,fb,fc;
  double min_x,min_y;

  a = 0.0001; 
  b = 1;
  _solve_erfc_alpha = alpha;
  ks_nr_mnbrake(&a,&b,&c,&fa,&fb,&fc,_solve_erfc);

  min_y = ks_nr_golden(a,b,c,_solve_erfc,1e-5,&min_x);
  return fabs(min_x-_solve_erfc_cut_off);
}
double ks_solve_erfc_alpha(double cut_off, double r)
{
  double a,b,c,fa,fb,fc;
  double min_x,min_y;

  _solve_erfc_sub = r;
  _solve_erfc_cut_off = cut_off;

  a = 0.001;
  b = 0.002;
  ks_nr_mnbrake(&a,&b,&c,&fa,&fb,&fc,_solve_erfc_alpha_func);

  min_y = ks_nr_golden(a,b,c,_solve_erfc_alpha_func,1e-5,&min_x);

  return min_x;
}
static double _solve_ewald_alpha;
static double _solve_ewald_ialpha2;
static double _solve_ewald_side;
static double _solve_ewald_iside2;
static double _solve_ewald_sub;
static double _solve_ewald_ivol;
static double _solve_ewald_kcut;
static double _solve_ewald_exp(double k)
{
  return fabs(exp(-KS_PI_PI*_solve_ewald_ialpha2*_solve_ewald_iside2*k*k)
	      *_solve_ewald_ivol/(k*k*_solve_ewald_iside2)
	      -_solve_ewald_sub);
}
double ks_solve_ewald_exp(double alpha, double side, double vol, double r)
{
  double a,b,c,fa,fb,fc;
  double min_x,min_y;

  _solve_ewald_sub = r;
  _solve_ewald_ivol = 1.0/vol;
  _solve_ewald_alpha = alpha;
  _solve_ewald_side = side;
  _solve_ewald_ialpha2 = 1.0/(_solve_ewald_alpha*_solve_ewald_alpha);
  _solve_ewald_iside2 = 1.0/(_solve_ewald_side*_solve_ewald_side);
  a = 0.1;
  b = 0.2;
  ks_nr_mnbrake(&a,&b,&c,&fa,&fb,&fc,_solve_ewald_exp);
  /*  printf("%f %f %f %f %f %f\n",a,b,c,fa,fb,fc);*/
  min_y = ks_nr_golden(a,b,c,_solve_ewald_exp,1e-5,&min_x);

  return min_x;
}
static double _solve_ewald_exp_alpha_func(double alpha)
{
  double a,b,c,fa,fb,fc;
  double min_x,min_y;

  _solve_ewald_alpha = alpha;
  _solve_ewald_ialpha2 = 1.0/(_solve_ewald_alpha*_solve_ewald_alpha);
  a = 0.1;
  b = 0.2;
  ks_nr_mnbrake(&a,&b,&c,&fa,&fb,&fc,_solve_ewald_exp);
  min_y = ks_nr_golden(a,b,c,_solve_ewald_exp,1e-5,&min_x);

  return fabs(min_x-_solve_ewald_kcut);
}
double ks_solve_ewald_exp_alpha(double k, double side, double vol, double r)
{
  double a,b,c,fa,fb,fc;
  double min_x,min_y;

  _solve_ewald_sub = r;
  _solve_ewald_ivol = 1.0/vol;
  _solve_ewald_side = side;
  _solve_ewald_iside2 = 1.0/(_solve_ewald_side*_solve_ewald_side);
  _solve_ewald_kcut = k;
  a = 0.1;
  b = 0.2;
  ks_nr_mnbrake(&a,&b,&c,&fa,&fb,&fc,_solve_ewald_exp_alpha_func);
  /*  printf("%f %f %f %f %f %f\n",a,b,c,fa,fb,fc);*/
  min_y = ks_nr_golden(a,b,c,_solve_ewald_exp_alpha_func,1e-5,&min_x);

  return min_x;
}

}

void EwaldBase::setSide(double side)
{
  // estimate alpha and k-cutoff from cutoff and side
  if (alpha == 0.0) {
    std::cout << "estimate alpha" << std::endl;
    alpha = ks_solve_erfc_alpha(cutoff,erfcCutValue);
    if (kCutoff == 0.0) {
      kCutoff = ks_solve_ewald_exp(alpha,side,side*side*side,ewaldExpCutValue);
    }
  }
  else if (kCutoff == 0.0) {
    kCutoff = side*cutoff*alpha*alpha/M_PI;
  }
  std::cout << "alpha " << alpha << "  kCutoff " << kCutoff << " " << side << std::endl;
//  exit(0);
}

namespace EwaldModule {
void estimate_alpha_kCutoff(double &alpha, double &kCutoff,
                            double side, double cutoff, 
                            double erfcCutValue,
                            double ewaldExpCutValue)
{
  if (alpha == 0.0) {
    alpha = ks_solve_erfc_alpha(cutoff,erfcCutValue);
    if (kCutoff == 0.0) {
      kCutoff = ks_solve_ewald_exp(alpha,side,side*side*side,ewaldExpCutValue);
    }
  }
  else if (kCutoff == 0.0) {
    kCutoff = side*cutoff*alpha*alpha/M_PI;
  }
}
double estimate_alpha(double cutoff, double erfcCutValue)
{
  return ks_solve_erfc_alpha(cutoff,erfcCutValue);
}
double estimate_kCutoff(double side, double cutoff, double alpha)
{
  return side*cutoff*alpha*alpha/M_PI;
}
}
