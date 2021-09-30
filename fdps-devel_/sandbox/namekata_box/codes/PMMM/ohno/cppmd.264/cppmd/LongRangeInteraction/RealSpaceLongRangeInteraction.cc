#include <config.h>

#include <cmath>
#include <iostream>
#include <limits>
#ifdef USE_FFTW2
# include <fftw.h>
# if 0
typedef fftw_complex  complex_t;
# else
struct complex_t : public fftw_complex {
  complex_t(): fftw_complex() {}
  complex_t(double real, double imag) : fftw_complex() {
    re = real;
    im = imag;
  }
  complex_t& operator*=(const double& x) {
    re *= x;
    im *= x;
    return *this;
  }
};
# endif
#else
# include <complex>
typedef std::complex<double>  complex_t;
# define c_re(c)  ((c).real())
# define c_im(c)  ((c).imag())
#endif
#include "RealSpaceLongRangeInteraction.h"

/*--------------------------------------------------------------------------*/

namespace {
const double sqrt_pi = sqrt(M_PI);
inline double sq(double x) { return x * x; }
inline double cube(double x) { return x * x * x; }
inline int center(int i, int number) {
  return (i > number / 2) ? i - number : i;
}
inline int xyzToIndex(const SpaceVector<int>& number_of_grid_point,
                      int x, int y, int z) {
  return x * (number_of_grid_point.y * number_of_grid_point.z)
      + y * (number_of_grid_point.z) + z;
}

template<typename T>
const SpaceVector<T> MulE(const SpaceVector<T> &a, const SpaceVector<T> &b) {
  return SpaceVector<T>(a.x * b.x, a.y * b.y, a.z * b.z);
}
template<typename T>
const SpaceVector<T> DivE(const SpaceVector<T> &a, const SpaceVector<T> &b) {
  return SpaceVector<T>(a.x / b.x, a.y / b.y, a.z / b.z);
}

// !!! these should be given as parameters
const double beta = 0.25;
const SpaceVector<int> number_of_grid_point(32, 32, 32);  // member of LongRangeInteraction
const Position grid_box_position(0.0, 0.0, 0.0);
const SpaceVector<double> grid_box_size(62.0, 62.0, 62.0);
// !!!
}  // namespace

/*--------------------------------------------------------------------------*/

RealSpaceChargeAssign::RealSpaceChargeAssign(int unitid,
                                             const LongRangeParameter& _param)
    : ChargeAssign(unitid, _param) {
}

void RealSpaceChargeAssign::assign(ParticleArray& particlearray,
                                   const std::vector<ParticleRange>& selfrange,
                                   ParticleArray& ghost,
                                   const std::vector<ParticleRange>& ghostrange,
                                   GridData& gridcharge)
{
  std::cout << "  ChargeAssign::assign " << unit_identifier << std::endl;

  const int xyz =
      number_of_grid_point.x * number_of_grid_point.y * number_of_grid_point.z;

  const SpaceVector<double> dx(grid_box_size.x / number_of_grid_point.x,
                               grid_box_size.y / number_of_grid_point.y,
                               grid_box_size.z / number_of_grid_point.z);
  const double beta2 = sq(beta);

  const double coeff = cube(beta / sqrt_pi);
  double charge_neutrality = 0.0;
  Position grid_point_position;
  for (int x = 0; x < number_of_grid_point.x; ++x) {
    grid_point_position.x = grid_box_position.x + x * dx.x;
    for (int y = 0; y < number_of_grid_point.y; ++y) {
      grid_point_position.y = grid_box_position.y + y * dx.y;
      for (int z = 0; z < number_of_grid_point.z; ++z) {
        grid_point_position.z = grid_box_position.z + z * dx.z;
        const int j = xyzToIndex(number_of_grid_point, x, y, z);
        gridcharge.gridvalue[j] = 0.0;
        std::vector<ParticleRange>::size_type it;
        for (it = 0; it < selfrange.size(); ++it) {     // self
          for (int i = selfrange[it].begin; i < selfrange[it].end; ++i) {
            const Particle& p = particlearray[i];
            const double r2 = (grid_point_position - p.position).norm2();
            gridcharge.gridvalue[j] += coeff * p.charge * exp(-r2 * beta2);
          }
        }
        for (it = 0; it < ghostrange.size(); ++it) {   // ghost
          for (int i = ghostrange[it].begin; i < ghostrange[it].end; ++i) {
            const Particle& p = ghost[i];
            const double r2 = (grid_point_position - p.position).norm2();
            gridcharge.gridvalue[j] += coeff * p.charge * exp(-r2 * beta2);
          }
        }
        charge_neutrality += gridcharge.gridvalue[j];
      }
    }
  }
  charge_neutrality /= xyz;
  for (int i = 0; i < xyz; ++i) {
    gridcharge.gridvalue[i] -= charge_neutrality;
  }
}

void RealSpaceChargeAssign::backinterpolate(
    ParticleArray& particlearray,
    const std::vector<ParticleRange>& selfrange,
    ParticleArray& ghost,
    const std::vector<ParticleRange>& ghostrange,
    GridData& gridpotential) {
  std::cout << "  ChargeAssign::backinterpolate " << unit_identifier
            << std::endl;

  const SpaceVector<double> dx(grid_box_size.x / number_of_grid_point.x,
                               grid_box_size.y / number_of_grid_point.y,
                               grid_box_size.z / number_of_grid_point.z);
  const double beta2 = sq(beta);

#if 0
  // direct part
  {
    std::vector<ParticleRange>::size_type it;
    for (it = 0; it < selfrange.size(); ++it) {
      for (int i = selfrange[it].begin; i < selfrange[it].end; ++i) {
        Particle& pi = particlearray[i];
        std::vector<ParticleRange>::size_type jt;
        for (jt = 0; jt < selfrange.size(); ++jt) {
          for (int j = selfrange[jt].begin; j < selfrange[jt].end; ++j) {
            Particle& pj = particlearray[j];

            const SpaceVector<double> vec = pi.position - pj.position;
            const double r2 = vec.norm2();
            if (r2 > 0.0) {
              const double r = sqrt(r2);
              const double ar = beta * r;
              const double qq = pi.charge * pj.charge;
              const double coeff = qq / (r*r*r) *
                  (erfc(ar) + 2.0*ar*exp(-sq(ar))/sqrt_pi);
              pi.force += vec * coeff;
            }
          }
        }
      }
    }
  }
#endif

  const double dxyz = dx.x * dx.y * dx.z;
  const double coeff_q = 2.0 * cube(beta / sqrt_pi) * beta2 * dxyz;
  Position grid_point_position;
  for (int x = 0; x < number_of_grid_point.x; ++x) {
    grid_point_position.x = grid_box_position.x + x * dx.x;
    for (int y = 0; y < number_of_grid_point.y; ++y) {
      grid_point_position.y = grid_box_position.y + y * dx.y;
      for (int z = 0; z < number_of_grid_point.z; ++z) {
        grid_point_position.z = grid_box_position.z + z * dx.z;
        const int j = xyzToIndex(number_of_grid_point, x, y, z);
        std::vector<ParticleRange>::size_type it;
        for (it = 0; it < selfrange.size(); ++it) {     // self
          for (int i = selfrange[it].begin; i < selfrange[it].end; ++i) {
            Particle& p = particlearray[i];
            const SpaceVector<double> vec = grid_point_position - p.position;
            const double r2 = vec.norm2();
            p.force -= gridpotential.gridvalue[j] * coeff_q * p.charge
                * exp(-r2 * beta2) * vec;
          }
        }
      }
    }
  }
#if 1
  for (int i = 0; i < 6; ++i) {
    printf("%d: (% g,% g,% g)\n",i,
           particlearray[i].force.x,
           particlearray[i].force.y,
           particlearray[i].force.z);
  }
#endif
}

/*--------------------------------------------------------------------------*/

RealSpacePoissonSolver::RealSpacePoissonSolver(int unitid, const LongRangeParameter& _param)
    : PoissonSolver(unitid,_param) {
}

void RealSpacePoissonSolver::solvePoisson(GridData& gridcharge,
                                          GridData& gridpotential,
                                          double& energy) {
  std::cout << "  PoissonSolver " << unit_identifier << std::endl;

  const int xyz =
      number_of_grid_point.x * number_of_grid_point.y * number_of_grid_point.z;

  const SpaceVector<double> dx(grid_box_size.x / number_of_grid_point.x,
                               grid_box_size.y / number_of_grid_point.y,
                               grid_box_size.z / number_of_grid_point.z);
  const SpaceVector<double> dx2(sq(dx.x), sq(dx.y), sq(dx.z));

  SpaceVector<double> dk(2.0 * M_PI / number_of_grid_point.x,
                         2.0 * M_PI / number_of_grid_point.y,
                         2.0 * M_PI / number_of_grid_point.z);

#ifdef USE_FFTW2
  fftwnd_plan plan = fftw3d_create_plan(
      number_of_grid_point.x, number_of_grid_point.y, number_of_grid_point.z,
      FFTW_BACKWARD, FFTW_MEASURE | FFTW_IN_PLACE);
  if (plan == NULL)
    return;   // !!! exit()
#endif  // USE_FFTW2

  complex_t *q = new complex_t[xyz];   // row major order (C order)
  for (int i = 0; i < xyz; ++i) {
    q[i] = complex_t(gridcharge.gridvalue[i],0.0);
  }

#ifdef USE_FFTW2
  fftwnd_one(plan, q, NULL);      // compute
#endif  // USE_FFTW2

  SpaceVector<double> pcoeff;
  for (int x = 0; x < number_of_grid_point.x; ++x) {
    pcoeff.x = (1 - cos(dk.x * center(x, number_of_grid_point.x))) / dx2.x;
    for (int y = 0; y < number_of_grid_point.y; ++y) {
      pcoeff.y = (1 - cos(dk.y * center(y, number_of_grid_point.y))) / dx2.y;
      for (int z = 0; z < number_of_grid_point.z; ++z) {
        pcoeff.z = (1 - cos(dk.z * center(z, number_of_grid_point.z))) / dx2.z;
        const double pcoeff_xyz_2 = pcoeff.x + pcoeff.y + pcoeff.z;
        const int i = xyzToIndex(number_of_grid_point, x, y, z);
        if (fabs(pcoeff_xyz_2) < std::numeric_limits<double>::epsilon()) {
          q[i] = complex_t(0.0,0.0);
        } else {  // pcoeff_xyz_2 != 0.0
          const double factor = 4.0 * M_PI / 2.0 / pcoeff_xyz_2;
          q[i] *= factor;
        }
      }
    }
  }

#ifdef USE_FFTW2
  fftwnd_one(plan, q, NULL);      // compute
#endif  // USE_FFTW2

#if 0
  double charge2 = 0.0;                     // !!! for self energy
#endif
  energy = 0.0;
  //const double beta2 = sq(beta);
  const double dxyz = dx.x * dx.y * dx.z;
  for (int x = 0; x < number_of_grid_point.x; ++x) {
    for (int y = 0; y < number_of_grid_point.y; ++y) {
      for (int z = 0; z < number_of_grid_point.z; ++z) {
        const int j = xyzToIndex(number_of_grid_point, x, y, z);
        q[j] = complex_t(c_re(q[j])/xyz,c_im(q[j]));
        energy += gridcharge.gridvalue[j] * c_re(q[j]);
        gridpotential.gridvalue[j] = c_re(q[j]);
#if 0
        std::vector<ParticleRange>::size_type it;
        for (it = 0; it < selfrange.size(); ++it) {     // self
          for (int i = selfrange[it].begin; i < selfrange[it].end; ++i) {
            charge2 += sq(particlearray[i].charge);      // !!! for self energy
          }
        }
        for (it = 0; it < ghostrange.size(); ++it) {   // ghost
          for (int i = ghostrange[it].begin; i < ghostrange[it].end; ++i) {
            charge2 += sq(ghost[i].charge);      // !!! for self energy
          }
        }
#endif
      }
    }
  }
  energy *= (0.5 * dxyz);

#if 0
  //energy += direct_part_potential;
  energy -= charge2 * beta / sqrt_pi;  // !!! for self energy
#endif

  delete []q;
#ifdef USE_FFTW2
  fftwnd_destroy_plan(plan);
#endif  // USE_FFTW2
}
