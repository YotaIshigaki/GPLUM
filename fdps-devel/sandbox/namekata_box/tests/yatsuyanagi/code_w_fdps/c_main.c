/* Standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
/* FDPS headers */
#include "user_defined.h"
#include "FDPS_c_if.h"
#include "parameter.h"

int setup_IC(int psys_num, int mark[MARK], int *mk)
{
    //   Get # of MPI processes and rank number
    int nprocs = fdps_get_num_procs();
    int myrank = fdps_get_rank();
    int i,j,cnt,num,rmn;
    double r,dr,theta;
    Full_particle *p;
    static double posX[MN], posY[MN];
    static double wz[MN];

    printf("myrank = %d\n", myrank);
    // Make an initial condition at RANK 0
    if (myrank == 0) {
        // ring
        dr = R0 / (double)MESH;
        cnt = 0;
        for (i=1; i<=MESH; i++) {
            r = i * dr;
            if (R1 <= r && r <= R2) {
                num = (int)floor(2.0 * PI * r / dr + 0.5);
                for (j=0; j<num; j++) {
                    theta = 2.0 * PI * j / ((double)num);
                    posX[cnt] = r * cos(theta);
                    posY[cnt] = r * sin(theta);
                    wz[cnt] = W0;
                    cnt++;
                    if (cnt >= MN) {
                        printf("Not enough array elements.");
                        return -1;
                    }
                }
            }
        }
        rmn = cnt;
        mark[*mk] = cnt;
        (*mk)++;
        printf("RING : %d\n", cnt);

        //Set # of local particles
        fdps_set_nptcl_loc(psys_num, rmn);
        Full_particle *p = (Full_particle *)fdps_get_psys_cptr(psys_num);

        for (i=0; i<rmn; i++) {
            p[i].pos.x = posX[i];
            p[i].pos.y = posY[i];
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION  // 3D mode is selected
            p[i].pos.z = 1.0;
#endif
            p[i].wz = wz[i];
        }
    } else {
        fdps_set_nptcl_loc(psys_num, 0);
    }

    return rmn;
}

// non mpi
int write_data(int psys_num, int t, int *mark)
{
    char fname[128];
    FILE *fp;
    int mk = 0;
    int i;
    int myrank = fdps_get_rank();
    Full_particle *p = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int lrmn = fdps_get_nptcl_loc(psys_num);

    sprintf(fname, "%03dw%02d_%02d.dat", t, mk, myrank);
    if ((fp = fopen(fname, "w")) == NULL) {
        printf("Cannot open %s\n", fname);
        return 1;
    }
    for (i=0; i<lrmn; i++) {
        if (i == mark[mk]) {
            fclose(fp);
            mk++;
            sprintf(fname, "%03dw%02d_%02d.dat", t, mk, myrank);
            if ((fp = fopen(fname, "w")) == NULL) {
                printf("Cannot open %s\n", fname);
                return 1;
            }
        }
        fprintf(fp, "%.10e %.10e %.4e\n", p[i].pos.x, p[i].pos.y, p[i].wz);
    }
    fclose(fp);

    return 0;
}

void check_data(int psys_num) {
    char fname[128];
    FILE *fp;
    sprintf(fname, "check.txt");
    if ((fp = fopen(fname,"w")) == NULL) {
        printf("Cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }
    int n_ptcl = fdps_get_nptcl_loc(psys_num);
    Full_particle *p = (Full_particle *) fdps_get_psys_cptr(psys_num);
    int i;
    for (i=0; i<n_ptcl; i++) {
        fprintf(fp, "%.15e   %.15e   %.15e   %.15e\n",
                p[i].pos.x, p[i].pos.y, p[i].vel.x, p[i].vel.y);
    }
    fclose(fp);
    exit(EXIT_SUCCESS);
}

int c_main()
{
    int Rmn;
    int error = 0;
    int t0, t1;
    int i;
    int mark[MARK], mk=0;
    Full_particle *p;
    static double posX[MN], posY[MN];
    static double u1X[MN], u1Y[MN], u2X[MN], u2Y[MN], u3X[MN], u3Y[MN], u4X[MN], u4Y[MN];

    printf("FDPS on C test code\n");

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
    printf("2D simlation mode is selected.\n");
#else
    printf("3D simlation mode is selected.\n");
#endif

    fdps_initialize();
    // Create and initialize dinfo object
    int dinfo_num;
    float coef_ema = 0.3;
    fdps_create_dinfo(&dinfo_num);
    fdps_init_dinfo(dinfo_num, coef_ema);
    // Create and initialize psys object
    int psys_num;
    fdps_create_psys(&psys_num, "full_particle");
    fdps_init_psys(psys_num);
    // Create and initialize tree object
    int tree_num;
    fdps_create_tree(&tree_num,
                     "Long,full_particle,full_particle,full_particle,Monopole");
    double theta = 0.5;
    //double theta = 0.0;
    //double theta = 0.001;
    int n_leaf_limit = 8;
    int n_group_limit = 64;
    fdps_init_tree(tree_num, MN, theta, n_leaf_limit, n_group_limit);
    // Make an initial condition
    mk = 0;
    for (i=0; i<MARK;i++) {
        mark[i] = -1;
    }
    if ((Rmn = setup_IC(psys_num, mark, &mk)) < 0) {
        error = 1;
        goto FINALIZE;
    }
    write_data(psys_num, 0, mark);
    // Domain decomposition and exchange particle
    fdps_decompose_domain_all(dinfo_num, psys_num, -1.0);
    fdps_exchange_particle(psys_num, dinfo_num);

    for (t0 = 1; t0 <= T0; t0++) {
        int lrmn = fdps_get_nptcl_loc(psys_num);
        printf("t0=%d lrmn=%d\n", t0,lrmn);
        fflush(stdout);
        for (t1=1; t1<= T1; t1++) {
            // 1st
            //printf(" CPU 1"); fflush(stdout);
            p = (Full_particle *)fdps_get_psys_cptr(psys_num);
            for (i=0; i<lrmn; i++) {
                posX[i] = p[i].pos.x;
                posY[i] = p[i].pos.y;
            }
//            cpu_Bs(cpu_posX, cpu_posY, cpu_u1X, cpu_u1Y, cpu_wz, rmn);
            fdps_calc_force_all_and_write_back(tree_num, calc_BS_ep_ep, calc_BS_ep_sp, psys_num, dinfo_num, true, FDPS_MAKE_LIST);
            //check_data(psys_num);
            for (i=0; i<lrmn; i++) {
                p[i].pos.x = posX[i] + 0.5 * DT * p[i].vel.x;
                p[i].pos.y = posY[i] + 0.5 * DT * p[i].vel.y;
                u1X[i] = p[i].vel.x;
                u1Y[i] = p[i].vel.y;
            }
            // 2nd
            //printf("  2"); fflush(stdout);
//            cpu_Bs(cpu_pos2X, cpu_pos2Y, cpu_u2X, cpu_u2Y, cpu_wz, rmn);
            fdps_calc_force_all_and_write_back(tree_num, calc_BS_ep_ep, calc_BS_ep_sp, psys_num, dinfo_num, true, FDPS_MAKE_LIST);
            for (i=0; i<lrmn; i++) {
                p[i].pos.x = posX[i] + 0.5 * DT * p[i].vel.x;
                p[i].pos.y = posY[i] + 0.5 * DT * p[i].vel.y;
                u2X[i] = p[i].vel.x;
                u2Y[i] = p[i].vel.y;
            }
            // 3rd
            //printf("  3"); fflush(stdout);
//            cpu_Bs(cpu_pos3X, cpu_pos3Y, cpu_u3X, cpu_u3Y, cpu_wz, rmn);
            fdps_calc_force_all_and_write_back(tree_num, calc_BS_ep_ep, calc_BS_ep_sp, psys_num, dinfo_num, true, FDPS_MAKE_LIST);
            for (i=0; i<lrmn; i++) {
                p[i].pos.x = posX[i] + 0.5 * DT * p[i].vel.x;
                p[i].pos.y = posY[i] + 0.5 * DT * p[i].vel.y;
                u3X[i] = p[i].vel.x;
                u3Y[i] = p[i].vel.y;
            }
            // 4th
            //printf("  4"); fflush(stdout);
//            cpu_Bs(cpu_pos4X, cpu_pos4Y, cpu_u4X, cpu_u4Y, cpu_wz, rmn);
            fdps_calc_force_all_and_write_back(tree_num, calc_BS_ep_ep, calc_BS_ep_sp, psys_num, dinfo_num, true, FDPS_MAKE_LIST);
            for (i=0; i<lrmn; i++) {
                u4X[i] = p[i].vel.x;
                u4Y[i] = p[i].vel.y;
                p[i].pos.x = posX[i] + DT * (u1X[i] + 2.0*u2X[i] + 2.0*u3X[i] + u4X[i]) / 6.0;
                p[i].pos.y = posY[i] + DT * (u1Y[i] + 2.0*u2Y[i] + 2.0*u3Y[i] + u4Y[i]) / 6.0;
            }
            check_data(psys_num);
        }
        fdps_decompose_domain_all(dinfo_num, psys_num, -1.0);
        fdps_exchange_particle(psys_num, dinfo_num);
        if (write_data(psys_num, t0, mark) != 0) {
            printf("failed to write data.\n");
            error = 1;
            goto FINALIZE;
        }
    }
    printf("calculation end\n");

FINALIZE:
    fdps_finalize();

    return error;
}

