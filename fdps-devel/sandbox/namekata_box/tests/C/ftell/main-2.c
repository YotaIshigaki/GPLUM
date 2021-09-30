#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct full_particle {
    long long id;
    double pos[3];
} Full_particle;

int main(int argc, char *argv[]) {
    // Local variables
    int i,j,k;

    // Set particle data to be output
    Full_particle ptcl_out;
    ptcl_out.id = 0;
    for (i=0; i<3; i++) {
       ptcl_out.pos[i] = (double)i;
    }
    fprintf(stdout,"sizeof(ptcl) = %lu [bytes]\n",sizeof(ptcl_out));

    // Output particle data
    FILE *fp;
    char filename[64];
    strcpy(filename,"data.dat");
    if ((fp = fopen(filename,"w")) == NULL) {
        fprintf(stderr,"cannot open file %s.\n",filename);
        exit(EXIT_FAILURE);
    }
    fprintf(fp,"%lld%lf%lf%lf\n",
            ptcl_out.id,ptcl_out.pos[0],ptcl_out.pos[1],ptcl_out.pos[2]);
    fclose(fp);

    // Read particle data and perform ftell
    Full_particle ptcl_in;
    if ((fp = fopen(filename,"r")) == NULL) {
        fprintf(stderr,"cannot open file %s.\n",filename);
        exit(EXIT_FAILURE);
    }
    fscanf(fp,"%lld%lf%lf%lf",
           &ptcl_in.id,&ptcl_in.pos[0],&ptcl_in.pos[1],&ptcl_in.pos[2]);
    long size_ptcl = ftell(fp);
    fprintf(stdout,"size_ptcl = %ld\n",size_ptcl);
    fclose(fp);

    return 0;
}
