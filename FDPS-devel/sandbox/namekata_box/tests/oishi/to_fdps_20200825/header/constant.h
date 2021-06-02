#pragma once


const PS::F64 DIM      = 3.0;
const PS::F64 dim_1    = 1.0 / DIM;
const PS::F64 pi       = atan(1.0)*4.0;
const PS::F64 SMTHK    = 1.2;//*dx
const PS::F64 dx       = 0.02;//initial particle spacing
const PS::F64 dx2      = dx * dx;
const PS::F64 kappa    = 2.0;//cubic spline...kappa=2.0
const PS::F64 smth_len = SMTHK * dx;// smoothing length = 1.2*d
const PS::F64 zero     = 0.0;
const PS::F64 one      = 1.0;
const PS::F64 two      = 2.0;
const PS::F64 gravity  = 9.80665;

/*Parameters*/
const PS::F64 dt              =   0.0001;
const PS::F64 end_time        =   1.0;
const PS::U32 OUTPUT_INTERVAL = 100;//100
//const PS::F64 C_MAX=0.1;
const PS::F64 rho0            = 1850.0;//initial density(kg/m3)
const PS::F64 rho_k           = 20.0 * 1000.0 / gravity;
const PS::F64 Height0         =  1.2;//H(m),
const PS::F64 Width0          = 19.0;//width(m)
const PS::F64 Depth0          =  7.5;
const PS::F64 Width1          =  4.0;
const PS::F64 Height1         =  5.0;
const PS::F64 Height2         =  6.2;
const PS::F64 Width2          = 10.0;

/***** itype list *****/
const PS::U32 S_W   = 10;//0-10:soil etc, 11-:wall(dummy)
const PS::U32 SOIL  =  0;//  SOIL PARTICLE
const PS::U32 KIBAN =  1;
// const PS::U32 UNDER = 11;//wall particle, under
// const PS::U32 LEFT  = 12;//wall particle, left
// const PS::U32 RIGHT = 13;
// const PS::U32 FRONT = 14;
// const PS::U32 BACK  = 15;
const PS::U32 WALL  = 20;//wall particle
// const PS::U32 GHOST =100;//
/*****************/

const PS::S32 WallLayer = 2;//number of wall layers(2+3)
const PS::F64 WallW     = Width0;//INNER wall width(boundary, not center of particle)
const PS::F64 WallH     = Height0 + 2.0 * dx;//INNER wall height
const PS::F64 WallH2    = Height2 + 2.0 * dx;
const PS::F64 WallD     = Depth0;
const PS::F64 Mass0     = rho0 * dx * dx * dx;//(kg)
const PS::F64 Mass_k    = rho_k * dx * dx * dx;
/* Parameters for kernel func. */
const PS::F64 val_W_2d     = 5.0 / (14.0 * pi * smth_len * smth_len);// 5/(14pi*h^2)
const PS::F64 val_W_3d     = 1.0 / (4.0 * pi * smth_len * smth_len * smth_len);// 1/(4*pi*h^3)
const PS::F64 val_gradW_2d = 15.0 / (14.0 * pi * smth_len * smth_len * smth_len * smth_len);// 15/(14pi*h^4)
const PS::F64 val_gradW_3d = 3.0 / (4.0 * pi * pow(smth_len, 5.0) );
const PS::F64 s_deld       = dx / smth_len;
/* Parameters for artificial viscosity*/
const PS::F64 alpha      = 0.01;
const PS::F64 beta       = 0.01;
const PS::F64 cvel       = 600.0;//sound speed in soil (m/s)
const PS::F64 multi_acij = - alpha * cvel;
/* Parameters for artificial stress */
const PS::F64 eps_AS     = 0.5;
const PS::F64 n_AS       = 2.55;
/* Parameters for motion eq. */


/*Parameters for stress rate*/
const PS::F64 E         = 1800.0*1000.0;//Young's elastic modulus(N/m2)
const PS::F64 E_k       = 8000.0*1000.0;//Young's elastic modulus(N/m2)
const PS::F64 Poiss     = 0.3;//Poisson's ratio, nu
// const PS::F64 K = E / (3.0 * (1.0 - 2.0*Poiss));//(N/m2=Pa)
// const PS::F64 G = 0.5 * E / (1.0 + Poiss);
// const PS::F64 Hook1 = (Poiss * E) / ((1.0 + Poiss) * (1.0 - 2.0 * Poiss));
// const PS::F64 Hook2 = E / (1.0 + Poiss);
// const PS::F64 EComp1 = 0.5 / G;
// const PS::F64 EComp2 = (1.0 - 2.0 * Poiss) / (3.0 * E);
const PS::F64 K         = 1.0 / (3.0 * (1.0 - 2.0*Poiss));// *Ey
const PS::F64 G         = 0.5 / (1.0 + Poiss);// *Ey
const PS::F64 Hook1     = Poiss / ((1.0 + Poiss) * (1.0 - 2.0 * Poiss));// *Ey
const PS::F64 Hook2     = 1.0 / (1.0 + Poiss);// *Ey
const PS::F64 EComp1    = 0.5 / G;// /Ey
const PS::F64 EComp2    = (1.0 - 2.0 * Poiss) / 3.0;// /Ey
//const PS::F64 coh0 = 9.8*1000.0;
const PS::F64 c1        = 5.0*1000.0;// Pa
const PS::F64 c_k       = 70.0*1000.0;
const PS::F64 cohR      = 0.0;
const PS::F64 phi       = 25.0;// internal friction [degree]
const PS::F64 phi_k     = 5.0;
const PS::F64 psi       = 0.0;// dilatancy angle
const PS::F64 alpha_c   = 0.0;//?
const PS::F64 tan_phi   = tan(phi * pi/180.0);
const PS::F64 sin_phi   = sin(phi * pi/180.0);
const PS::F64 cos_phi   = cos(phi * pi/180.0);
const PS::F64 sin_phi_k = sin(phi_k * pi/180.0);
const PS::F64 cos_phi_k = cos(phi_k * pi/180.0);
const PS::F64 tan_phi_k = tan(phi_k * pi/180.0);
const PS::F64 sin_psi   = sin(psi * pi/180.0);
//const PS::F64 al_phi = tan_phi/(sqrt(9.0+12.0*tan_phi*tan_phi));
//const PS::F64 kc = 3.0/(sqrt(9.0+12.0*tan_phi*tan_phi));
// const PS::F64mat del_2d(1.0, 1.0, 0.0);//Kronecker delta 2D
const PS::F64mat del_3d(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
const PS::S32 ELASTIC = 2;
const PS::S32 PLASTIC = 3;

/************ output file for C++ *****************/
std::ofstream ene("energy.csv");


