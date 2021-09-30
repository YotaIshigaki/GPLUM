#include <iomanip>
#include "Common.h"

#if 0
#define _GNU_SOURCE 1
#include <fenv.h>

static void __attribute__ ((constructor))
trapfpe ()
{
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#endif

double calcgrid(const SpaceVector<double>& box,
                const SpaceVector<int>& gridnums,
                const double alpha,
                const std::vector<Position>& cd,
                const std::vector<double>& charge,
                const int nexternal=1)
{
  std::cout << "box  " << box << std::endl;
  std::cout << "grid " << gridnums << std::endl;
  double mgalpha = sqrt(2.0)*alpha;
  std::cout << "alpha   " << alpha << std::endl;
  std::cout << "mgalpha " << mgalpha << std::endl;
  std::vector<Position> shift;
  for (int ix=-nexternal;ix<=nexternal;++ix) {
    for (int iy=-nexternal;iy<=nexternal;++iy) {
      for (int iz=-nexternal;iz<=nexternal;++iz) {
        shift.push_back(Position(box.x*ix, box.y*iy, box.z*iz));
      }
    }
  }
  SpaceVector<double> h(box.x/gridnums.x,
                        box.y/gridnums.y,
                        box.z/gridnums.z);
  Position gpos;
  double erec = 0.0;
  double erecm = 0.0;
  for (int ix=0;ix<gridnums.x;++ix) {
    gpos.x = ix*h.x;
    std::cerr << ix << "\r";
    for (int iy=0;iy<gridnums.y;++iy) {
      gpos.y = iy*h.y;
      for (int iz=0;iz<gridnums.z;++iz) {
        gpos.z = iz*h.z;
        double gpot = 0.0;
        double gc = 0.0;
        for (std::vector<Position>::const_iterator it = shift.begin();
             it != shift.end();++it) {
          Position gspos(gpos-(*it));
          for (int i=0;i<cd.size();++i) {
            double r = (gspos-cd[i]).norm();
            gc += charge[i]*exp(-mgalpha*mgalpha*r*r);
            gpot += charge[i]*erf(mgalpha*r)/r;
          }
        }
        gc *= mgalpha*mgalpha*mgalpha/sqrt(M_PI*M_PI*M_PI);
        erec += gc*gpot;
#ifdef MG_DUMP
        if (std::abs(gc) > 0.2) {
          erecm += gc*gpot;
          std::cout << "M " << ix << " " << iy << " " << iz << " " << gc << " " << gpot << " " << erecm << std::endl;
        }
#endif
      }
    }
  }
  erec *= h.x*h.y*h.z*0.5;
  std::cout << "Erec  " << erec << std::endl;
#ifdef MG_DUMP
  erecm *= h.x*h.y*h.z*0.5;
  std::cout << "ErecM  " << erecm << std::endl;
#endif
  double realpot = 0.0;
  double self = 0.0;
  for (int i=0;i<cd.size();++i) {
    self += charge[i]*charge[i];
    for (std::vector<Position>::const_iterator it = shift.begin();
         it != shift.end();++it) {
      Position cds = cd[i]-(*it);
      for (int j=0;j<cd.size();++j) {
        double r = (cds-cd[j]).norm();
        if (r > 0.0) {
          realpot += charge[i]*charge[j]*erfc(alpha*r)/r;
        }
      }
    }
  }
  realpot *= 0.5;
  std::cout << "Ereal " << realpot << std::endl;
  self *= -alpha*M_2_SQRTPI*0.5;
  std::cout << "Eself " << self << std::endl;
  double gridpot = realpot+erec+self;
  return gridpot;
}

double calcdirect(const SpaceVector<double>& box,
                  const std::vector<Position>& cd,
                  const std::vector<double>& charge,
                  const int nexternal=1)
{
  std::vector<Position> shift;
  for (int ix=-nexternal;ix<=nexternal;++ix) {
    for (int iy=-nexternal;iy<=nexternal;++iy) {
      for (int iz=-nexternal;iz<=nexternal;++iz) {
        shift.push_back(Position(box.x*ix, box.y*iy, box.z*iz));
      }
    }
  }
  double directpot = 0;
  for (std::vector<Position>::const_iterator it = shift.begin();
       it != shift.end();++it) {
    for (int i=0;i<cd.size();++i) {
      Position pos = cd[i] + (*it);
      for (int j=0;j<cd.size();++j) {
        double r = (cd[j]-pos).norm();
        if (r > 0.0) {
          directpot += charge[i]*charge[j]/r;
        }
      }
    }
  }
  directpot *= 0.5;
  return directpot;
}
              
void test(int gridnum, int nexternal, double w)
{
  if (w == 0.0) {
    w = 50.0;
  }
  double tt = 2.0;
  //double alpha = sqrt(0.5);
  double alpha = 1.0;
  std::vector<Position> cd;
  cd.push_back(Position(0.5*w+0.5*tt,0.5*w,0.5*w));
  cd.push_back(Position(0.5*w-0.5*tt,0.5*w,0.5*w));
  std::vector<double> charge;
  charge.push_back(3.0);
  charge.push_back(-3.0);
  SpaceVector<double> box(w);
  SpaceVector<int> gridnums(gridnum);
  double gridpot = calcgrid(box, gridnums, alpha, cd, charge, nexternal);
  double directpot = calcdirect(box, cd, charge, nexternal);
  std::cout << std::setprecision(16);
  std::cout << "grid   = " << gridpot << std::endl;
  std::cout << "direct = " << directpot << std::endl;
}

void test2(int gridnum, int nexternal, int direction=0)
{
  SpaceVector<double> L(62.0, 62.0, 62.0);
  const double bond = 0.957;
  const double angle = 104.5;
  const double theta = M_PI / 180.0 * angle / 2.0;
  const double bond_sin_theta = bond * sin(theta);
  const double bond_cos_theta = bond * cos(theta);
  double alpha = 0.25;
  SpaceVector<int> gridnums(gridnum);

  std::vector<Position> cd(6);
  std::vector<double> charge(6);

  // water 1 -- Oxygen
  cd[0] = Position(L.x/2, L.y/2, L.z/2 + 1.0);
  charge[0] = -0.82;
  // water 1 -- Hydrogen
  cd[1] = cd[0];
  cd[1].x += bond_sin_theta;
  cd[1].z += bond_cos_theta;
  charge[1] = 0.41;
  // water 1 -- Hydrogen
  cd[2] = cd[0];
  cd[2].x -= bond_sin_theta;
  cd[2].z += bond_cos_theta;
  charge[2] = 0.41;

  // water 2 -- Oxygen
  cd[3] = Position(L.x/2, L.y/2, L.z/2 - 1.0);
  charge[3] = -0.82;
  // water 2 -- Hydrogen
  cd[4] = cd[3];
  if (direction == 0) {
    cd[4].x += bond_sin_theta;
    cd[4].z -= bond_cos_theta;
  }
  else {
    cd[4].x -= bond_sin_theta;
    cd[4].z += bond_cos_theta;
  }
  charge[4] = 0.41;
  // water 2 -- Hydrogen
  cd[5] = cd[3];
  if (direction == 0) {
    cd[5].x -= bond_sin_theta;
    cd[5].z -= bond_cos_theta;
  }
  else {
    cd[5].x += bond_sin_theta;
    cd[5].z += bond_cos_theta;
  }
  charge[5] = 0.41;

  double gridpot = calcgrid(L, gridnums, alpha, cd, charge, nexternal);
  double directpot = calcdirect(L, cd, charge, nexternal);
  std::cout << std::setprecision(16);
  std::cout << "Egrid   = " << gridpot << std::endl;
  std::cout << "Edirect = " << directpot << std::endl;
}

int main(int argc, char *argv[])
{
  int gridnum = 128;
  int nexternal = 1;
  double w = 0.0;
  int direction = 0;
  if (argc > 1) {
    gridnum = atoi(argv[1]);
  }
  if (argc > 2) {
    nexternal = atoi(argv[2]);
  }
  if (argc > 3) {
    w = atof(argv[3]);
  }
  if (argc > 4) {
    direction = atoi(argv[4]);
  }
  test(gridnum, nexternal, w);
  //test2(gridnum, nexternal, direction);
}
