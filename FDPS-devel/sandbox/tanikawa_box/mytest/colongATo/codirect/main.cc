#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

typedef struct nbody {
    int    id;
    double mass;
    double pos[3];
    double vel[3];
    double acc[3];
    double pot;
    int njreal;
}NBODY;

void loadInitialCondition(FILE *fp,
                          int *n,
                          NBODY *system)
{
    int itmp;
    double dtmp;

    fscanf(fp, "%lf", &dtmp);
    fscanf(fp, "%d", n);
    
    for(int i = 0; i < *n; i++) {
        system[i].id     = i;
        fscanf(fp, "%d", &itmp);
        fscanf(fp, "%lf", &(system[i].mass));
        fscanf(fp, "%lf%lf%lf", &(system[i].pos[0]), &(system[i].pos[1]), &(system[i].pos[2]));
        fscanf(fp, "%lf%lf%lf", &(system[i].vel[0]), &(system[i].vel[1]), &(system[i].vel[2]));
    }

    return;
}

int main(int argc, char **argv)
{
    double eps2 = 1. / (32. * 32.);
    double cutoff2 = 1. / (8. * 8.) + eps2;
    int n;
    NBODY system[32768];

    FILE *fp = fopen(argv[1], "r");
    loadInitialCondition(fp, &n, system);
    fclose(fp);

    for(int i = 0; i < n; i++) {
        double posi[3] = {system[i].pos[0], system[i].pos[1], system[i].pos[2]};
        double acci[3] = {0., 0., 0.};
        double poti    = 0.;
        int njreal =0;
        for(int j = 0; j < n; j++) {
            double mj      = system[j].mass;
            double posj[3] = {system[j].pos[0], system[j].pos[1], system[j].pos[2]};
            double dx[3] = {posj[0] - posi[0], posj[1] - posi[1], posj[2] - posi[2]};

            double r2   = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] + eps2;
            double rinv = ((j != i) ? 1. / sqrt(r2) : 0.);
            rinv        = ((r2 < cutoff2) ? rinv : 0.);
            njreal     += ((r2 < cutoff2) ? 1 : 0);
            double pot  = mj * rinv;
            poti -= pot;

            double fij = pot * rinv * rinv;
            acci[0] += fij * dx[0];
            acci[1] += fij * dx[1];
            acci[2] += fij * dx[2];
        }
        system[i].acc[0] = acci[0];
        system[i].acc[1] = acci[1];
        system[i].acc[2] = acci[2];
        system[i].pot    = poti;
        system[i].njreal = njreal;
    }

    {
        FILE *fp = fopen("dump.dat", "w");
        for(int i = 0; i < n; i++) {
            fprintf(fp, "%6d", system[i].id);
            fprintf(fp, " %+e %+e %+e",
                    system[i].acc[0], system[i].acc[1], system[i].acc[2]);
            fprintf(fp, " %+e", system[i].pot);
            fprintf(fp, " %6d\n", system[i].njreal);
        }
        fclose(fp);
    }

    return 0;
}
