#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// strength of the vortex
#define W0 0.15

// particle density
#define MESH 25

// radius of circular wall
#define R0 50.0

// inner radius of ring
#define R1 3.0

// outer radius of ring
#define R2 40.0

#define PI 3.14159265358979

#define MN 100000 // max number of particles

#define DT 1.0e-4
#define T0 10
#define T1 1

#define SQR(x) ((x)*(x))

int initializer(double posX[], double posY[], double wz[], int cnt);
int init_ring(double posX[], double posY[], double wz[], int cnt);
void cpu_Bs2(double posX[], double posY[], double uX[], double uY[], double wz[MN], int rmn);
int write_data(double posX[], double posY[], double wz[], int mn, int t);

int initializer(double posX[], double posY[], double wz[], int cnt)
{
	return init_ring(posX, posY, wz, cnt);
}

int init_ring(double posX[], double posY[], double wz[], int cnt)
{
	int i,j;
	double r, dr, theta;
	double x, y;
	int num;

	printf("n: MESH = %d\n", MESH);
	printf("q: R2   = %f  (outer radi of ring)\n", R2);
	printf("v: R1   = %f  (radi of core)\n", R1);

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
	printf("RING : %d\n", cnt);

	return cnt;
}


int main(void)
{
	double *cpu_posX = NULL, *cpu_pos2X = NULL, *cpu_pos3X = NULL, *cpu_pos4X = NULL;
	double *cpu_posY = NULL, *cpu_pos2Y = NULL, *cpu_pos3Y = NULL, *cpu_pos4Y = NULL;
	double *cpu_u1X = NULL, *cpu_u2X = NULL, *cpu_u3X = NULL, *cpu_u4X = NULL;
	double *cpu_u1Y = NULL, *cpu_u2Y = NULL, *cpu_u3Y = NULL, *cpu_u4Y = NULL;
	double *cpu_wz = NULL;
	int cnt;
	int rmn;
   size_t mem_size;

	static double ini_posX[MN];
	static double ini_posY[MN];
	static double ini_wz[MN];
	int t0,t1;
	int i,j;

	cnt = 0;
	// initial profile
	if ((cnt = initializer(ini_posX, ini_posY, ini_wz, cnt)) < 0) {
		goto ERROR_ALL;
	}
	rmn = cnt;
	if (rmn > MN) {
		printf("Not enough memory for psudo particles. \n");
		goto ERROR_ALL;
	}

	//memory allocate
	mem_size = sizeof(double) * rmn;
	cpu_posX = (double *)malloc(mem_size);
	if (cpu_posX == NULL) {
		printf("malloc() for cpu_posX failed.\n");
		goto ERROR_ALL;
	}
	cpu_posY = (double *)malloc(mem_size);
	if (cpu_posY == NULL) {
		printf("malloc() for cpu_posY failed.\n");
		goto ERROR_ALL;
	}
	cpu_pos2X = (double *)malloc(mem_size);
	if (cpu_pos2X == NULL) {
		printf("malloc() for cpu_posX failed.\n");
		goto ERROR_ALL;
	}
	cpu_pos2Y = (double *)malloc(mem_size);
	if (cpu_pos2Y == NULL) {
		printf("malloc() for cpu_posY failed.\n");
		goto ERROR_ALL;
	}
	cpu_pos3X = (double *)malloc(mem_size);
	if (cpu_pos3X == NULL) {
		printf("malloc() for cpu_posX failed.\n");
		goto ERROR_ALL;
	}
	cpu_pos3Y = (double *)malloc(mem_size);
	if (cpu_pos3Y == NULL) {
		printf("malloc() for cpu_posY failed.\n");
		goto ERROR_ALL;
	}
	cpu_pos4X = (double *)malloc(mem_size);
	if (cpu_pos4X == NULL) {
		printf("malloc() for cpu_posX failed.\n");
		goto ERROR_ALL;
	}
	cpu_pos4Y = (double *)malloc(mem_size);
	if (cpu_pos4Y == NULL) {
		printf("malloc() for cpu_posY failed.\n");
		goto ERROR_ALL;
	}
	cpu_u1X = (double *)malloc(mem_size);
	if (cpu_u1X == NULL) {
		printf("malloc() for cpu_uX failed.\n");
		goto ERROR_ALL;
	}
	cpu_u1Y = (double *)malloc(mem_size);
	if (cpu_u1Y == NULL) {
		printf("malloc() for cpu_uY failed.\n");
		goto ERROR_ALL;
	}
	cpu_u2X = (double *)malloc(mem_size);
	if (cpu_u2X == NULL) {
		printf("malloc() for cpu_uX failed.\n");
		goto ERROR_ALL;
	}
	cpu_u2Y = (double *)malloc(mem_size);
	if (cpu_u2Y == NULL) {
		printf("malloc() for cpu_uY failed.\n");
		goto ERROR_ALL;
	}
	cpu_u3X = (double *)malloc(mem_size);
	if (cpu_u3X == NULL) {
		printf("malloc() for cpu_uX failed.\n");
		goto ERROR_ALL;
	}
	cpu_u3Y = (double *)malloc(mem_size);
	if (cpu_u3Y == NULL) {
		printf("malloc() for cpu_uY failed.\n");
		goto ERROR_ALL;
	}
	cpu_u4X = (double *)malloc(mem_size);
	if (cpu_u4X == NULL) {
		printf("malloc() for cpu_uX failed.\n");
		goto ERROR_ALL;
	}
	cpu_u4Y = (double *)malloc(mem_size);
	if (cpu_u4Y == NULL) {
		printf("malloc() for cpu_uY failed.\n");
		goto ERROR_ALL;
	}
	cpu_wz = (double *)malloc(mem_size);
	if (cpu_wz == NULL) {
		printf("malloc() for cpu_wz failed.\n");
		goto ERROR_ALL;
	}
	for (i=0; i<rmn; i++) {
		cpu_posX[i] = ini_posX[i];
		cpu_posY[i] = ini_posY[i];
		cpu_wz[i] = ini_wz[i];
	}
	printf("malloc end.\n");

	if (write_data(cpu_posX, cpu_posY,cpu_wz, rmn, 0) != 0) {
		printf("failed to write data.\n");
		goto ERROR_ALL;
	}

	for (t0 = 1; t0 <= T0; t0++) {
		printf("t0=%d\n", t0);
		fflush(stdout);
		for (t1=1; t1<= T1; t1++) {
			// 1st
			cpu_Bs2(cpu_posX, cpu_posY, cpu_u1X, cpu_u1Y, cpu_wz, rmn);
			for (i=0; i<rmn; i++) {
				cpu_pos2X[i] = cpu_posX[i] + 0.5 * DT * cpu_u1X[i];
				cpu_pos2Y[i] = cpu_posY[i] + 0.5 * DT * cpu_u1Y[i];
			}
			// 2nd
			cpu_Bs2(cpu_pos2X, cpu_pos2Y, cpu_u2X, cpu_u2Y, cpu_wz, rmn);
			for (i=0; i<rmn; i++) {
				cpu_pos3X[i] = cpu_posX[i] + 0.5 * DT * cpu_u2X[i];
				cpu_pos3Y[i] = cpu_posY[i] + 0.5 * DT * cpu_u2Y[i];
			}
			// 3rd
			cpu_Bs2(cpu_pos3X, cpu_pos3Y, cpu_u3X, cpu_u3Y, cpu_wz, rmn);
			for (i=0; i<rmn; i++) {
				cpu_pos4X[i] = cpu_posX[i] + DT * cpu_u3X[i];
				cpu_pos4Y[i] = cpu_posY[i] + DT * cpu_u3Y[i];
			}
			// 4th
			cpu_Bs2(cpu_pos4X, cpu_pos4Y, cpu_u4X, cpu_u4Y, cpu_wz, rmn);
			for (i=0; i<rmn; i++) {
				cpu_posX[i] += DT * (cpu_u1X[i] + 2.0*cpu_u2X[i] + 2.0*cpu_u3X[i] + cpu_u4X[i]) / 6.0;
				cpu_posY[i] += DT * (cpu_u1Y[i] + 2.0*cpu_u2Y[i] + 2.0*cpu_u3Y[i] + cpu_u4Y[i]) / 6.0;
			}
            printf("%.15e %.15e\n",cpu_posX[0],cpu_posY[0]);
            exit(0);
		}
		if (write_data(cpu_posX, cpu_posY, cpu_wz, rmn, t0) != 0) {
			printf("failed to write data.\n");
			goto ERROR_ALL;
		}
	}
	printf("calculation end\n");

ERROR_ALL:
	// final handler
	free(cpu_posX);
	free(cpu_pos2X);
	free(cpu_pos3X);
	free(cpu_pos4X);
	free(cpu_posY);
	free(cpu_pos2Y);
	free(cpu_pos3Y);
	free(cpu_pos4Y);
	free(cpu_u1X);
	free(cpu_u2X);
	free(cpu_u3X);
	free(cpu_u4X);
	free(cpu_u1Y);
	free(cpu_u2Y);
	free(cpu_u3Y);
	free(cpu_u4Y);
	free(cpu_wz);

    return 0;
}

// Special thanks to Makino-san. ^^)
void cpu_Bs2(double posX[], double posY[], double uX[], double uY[], double wz[MN], int rmn)
{
	int i, j;
	double r2;

	for (i=0; i<rmn; i++) {
		uX[i] = 0.0;
		uY[i] = 0.0;
		for (j=0; j<i; j++) {
			// real
			const double posXij = posX[i] - posX[j];
			const double posYij = posY[i] - posY[j];
			const double rij2 = posXij * posXij + posYij * posYij;
			// pseudo
			const double rj2 = posX[j] * posX[j] + posY[j] * posY[j];
			const double pposXj = posX[j] * R0 * R0 / rj2;
			const double pposYj = posY[j] * R0 * R0 / rj2;
			const double pposXij = posX[i] - pposXj;
			const double pposYij = posY[i] - pposYj;
			const double prij2 = pposXij * pposXij + pposYij * pposYij;
			uX[i] -= wz[j] * posYij / rij2 - wz[j] * pposYij / prij2;
			uY[i] += wz[j] * posXij / rij2 - wz[j] * pposXij / prij2;
		}
		// i=j pseudo only
		const double rj2 = posX[j] * posX[j] + posY[j] * posY[j];
		const double pposXj = posX[j] * R0 * R0 / rj2;
		const double pposYj = posY[j] * R0 * R0 / rj2;
		const double pposXij = posX[i] - pposXj;
		const double pposYij = posY[i] - pposYj;
		const double prij2 = pposXij * pposXij + pposYij * pposYij;
		uX[i] -= - wz[j] * pposYij / prij2;
		uY[i] += - wz[j] * pposXij / prij2;
		//
		for (j=i+1; j<rmn; j++) {
			// real
			const double posXij = posX[i] - posX[j];
			const double posYij = posY[i] - posY[j];
			const double rij2 = posXij * posXij + posYij * posYij;
			// pseudo
			const double rj2 = posX[j] * posX[j] + posY[j] * posY[j];
			const double pposXj = posX[j] * R0 * R0 / rj2;
			const double pposYj = posY[j] * R0 * R0 / rj2;
			const double pposXij = posX[i] - pposXj;
			const double pposYij = posY[i] - pposYj;
			const double prij2 = pposXij * pposXij + pposYij * pposYij;
			uX[i] -= wz[j] * posYij / rij2 - wz[j] * pposYij / prij2;
			uY[i] += wz[j] * posXij / rij2 - wz[j] * pposXij / prij2;
		}
	}
}
int write_data(double posX[], double posY[], double wz[], int rmn, int t)
{
	char fname[128];
	FILE *fp;
	int i, mk = 0;
	double r2, d;

	sprintf(fname, "t%03d.dat", t);
	if ((fp = fopen(fname, "w")) == NULL) {
		printf("Cannot open %s\n", fname);
		return 1;
	}
	for (i=0; i<rmn; i++) {
		fprintf(fp, "%.15e %.15e\n", posX[i], posY[i]);
	}
	fclose(fp);

	return 0;
}

