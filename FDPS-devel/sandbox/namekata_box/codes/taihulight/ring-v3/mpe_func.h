#pragma once

struct EpiMM{
  double x,y,z,mass;
  double vx,vy,vz;
  long id;
};
 
#define EpjMM EpiMM

struct SpjMM{
  double x,y,z,mass;
};

struct ForceMM{
  double ax,ay,az,pot;
};

void calcForcePlanet(const EpiMM epi[], const int n_epi, const double mp, ForceMM force[]);

void calcDirectGrav(const EpiMM epi[], const int n_epi,
                    const EpjMM epj[], const int n_epj, const int adr_epj[],
                    ForceMM  force[], double r_epi, double r_epj, double kappa, double eta, const bool adr_flag=true);

void calcSPGrav(const EpiMM epi[], const int n_epi,
                const SpjMM spj[], const int n_spj, const int adr_spj[],
                ForceMM  force[], const bool adr_flag=true);

struct MPE_pars{
  double kappa,eta,r_ep,r_sat;
};
