#include <stdio.h>
#include <stdlib.h>
#include <string.h>

inline double frand() {
    return (double) rand() / ((double)RAND_MAX + 1.);
}

void makeForGreeM(FILE *fp);
void outputGreeM(char *inputfile);

int main(int argc,
         char **argv) {

    /*
    int n = 8192;

    printf("%.1f\n", 0.0);
    printf("%d\n", n);

    for(int i = 0; i < n; i++) {
        printf("%d\n", 1);
        printf("%+.16e\n", 1. / (double)n);
        double x = frand();
        double y = frand();
        double z = frand();
        printf("%+.16e %+.16e %+.16e\n", x, y, z);
        printf("%+.16e %+.16e %+.16e\n", 0., 0., 0.);
    }
    */

    if(argc != 3) {
        fprintf(stderr, "%s <file> <input/output>\n", argv[0]);
        exit(0);
    }
    if(strcmp(argv[2], "input") == 0) {
        FILE *fp = fopen(argv[1], "w");
        makeForGreeM(fp);
        fclose(fp);
    } else if(strcmp(argv[2], "output") == 0) {
        outputGreeM(argv[1]);
    } else {
        fprintf(stderr, "arg2 should be input/output.\n");
        exit(0);
    }
    
    return 0;
}

#define NDM (32768)
void makeForGreeM(FILE *fp)
{
    int n = NDM;
    int itmp    = 0;
    float ftmp  = 1.;
    double dtmp = 1.;
    
    fwrite(&itmp, sizeof(int),    1, fp);
    fwrite(&n,    sizeof(int),    1, fp);
    fwrite(&itmp, sizeof(int),    1, fp);
    ftmp = 0.3; // omega0
    fwrite(&ftmp, sizeof(float),  1, fp);
    ftmp = 0.1; // omageb
    fwrite(&ftmp, sizeof(float),  1, fp);
    ftmp = 0.7; // lambda0
    fwrite(&ftmp, sizeof(float),  1, fp);
    ftmp = 0.7; // hubble
    fwrite(&ftmp, sizeof(float),  1, fp);
    fwrite(&ftmp, sizeof(float),  1, fp);
    fwrite(&ftmp, sizeof(float),  1, fp);
    fwrite(&ftmp, sizeof(float),  1, fp);
    fwrite(&dtmp, sizeof(double), 1, fp);
    fwrite(&dtmp, sizeof(double), 1, fp);
    fwrite(&dtmp, sizeof(double), 1, fp);

    static float x[NDM][3], v[NDM][3];
    static long long int id[NDM];

    for(int i = 0; i < n; i++) {
        x[i][0] = frand();
        x[i][1] = frand();
        x[i][2] = frand();
        v[i][0] = 0.;
        v[i][1] = 0.;
        v[i][2] = 0.;
        id[i]   = i;
    }

    fwrite(x, sizeof(float), n*3, fp);
    fwrite(v, sizeof(float), n*3, fp);
    fwrite(id, sizeof(long long int), n, fp);
    
    return;
}

void outputGreeM(char *inputfile)
{
    static long long int id[NDM];
    static float x[NDM][3], v[NDM][3];

    int n;
    int itmp;
    float ftmp;
    double dtmp;

    FILE *fp = fopen(inputfile, "r");
    
    fread(&itmp, sizeof(int),    1, fp);
    fread(&n,    sizeof(int),    1, fp);
    fread(&itmp, sizeof(int),    1, fp);
    fread(&itmp, sizeof(float),  1, fp);
    fread(&itmp, sizeof(float),  1, fp);
    fread(&itmp, sizeof(float),  1, fp);
    fread(&itmp, sizeof(float),  1, fp);
    fread(&itmp, sizeof(float),  1, fp);
    fread(&itmp, sizeof(float),  1, fp);
    fread(&itmp, sizeof(float),  1, fp);
    fread(&itmp, sizeof(double), 1, fp);
    fread(&itmp, sizeof(double), 1, fp);
    fread(&itmp, sizeof(double), 1, fp);

    fread(x, sizeof(float), n*3, fp);
    fread(v, sizeof(float), n*3, fp);
    fread(id, sizeof(long long int), n, fp);

    fclose(fp);

    for(int i = 0; i < n; i++) {
        printf("%8d %+.10e %+.10e %+.10e", id[i], x[i][0], x[i][1], x[i][2]);
//        printf("%+e %+e %+e", x[i][0], x[i][1], x[i][2]);
//        printf(" %+e %+e %+e", v[i][0], v[i][1], v[i][2]);
        printf("\n");
    }
    
    return;
}

