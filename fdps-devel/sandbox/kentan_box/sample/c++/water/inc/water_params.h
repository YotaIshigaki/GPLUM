#ifndef H_WATER_PARAMS
#define H_WATER_PARAMS

#include "unit.h"

//SPC/E parameters
const double MASS_OXY    = 15.9994 * 1e-3 / Na / unit_mass;
const double SIGMA_OXY   = 3.166 * 1e-10 / unit_length;
const double EPSILON_OXY = 0.6500 * 1e+3 / Na / unit_energy;
const double CHARGE_OXY  = -0.8476 * ele / unit_coulomb;

const double MASS_HYD    = 1.00800 * 1e-3 / Na / unit_mass;
const double SIGMA_HYD   = 0.0;
const double EPSILON_HYD = 0.0;
const double CHARGE_HYD  = 0.4238 * ele / unit_coulomb;

const double ANGLE_HOH   = M_PI / 180.0 * 109.47;

const double BOND_OH = 1.0 * 1e-10 / unit_length;
const double BOND_HH = 2.0 * BOND_OH * cos(M_PI/180.0*(90.0 - 0.5*109.47));

const double BOND_COEF  = 3450.;
const double ANGLE_COEF = 383.;
#endif
