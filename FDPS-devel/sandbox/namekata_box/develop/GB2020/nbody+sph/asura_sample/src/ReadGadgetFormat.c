#include "config.h"

static void do_what_you_want(void);
static void Gadget_unit_conversion(void);
static void Gadget_allocate_memory(void);
static void Gadget_load_snapshot(char *fname, int files);
static void Gadget_reordering(void);

struct Gadget_io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} GadgetHeader1;

int     NumPart, Ngas;

struct Gadget_particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} *GadgetP;

int *Id;

double  Time, Redshift;

/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int ReadGadetFormat(void){

  char path[MaxCharactersInLine], input_fname[MaxCharactersInLine], basename[MaxCharactersInLine];

  sprintf(path, "./");
  sprintf(basename, "data");
  int snapshot_number= 6;                    /* number of snapshot */
  int files=1;                               /* number of files per snapshot */

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  Gadget_load_snapshot(input_fname, files);

  Gadget_reordering();  /* call this routine only if your ID's are set properly */

  Gadget_unit_conversion();  /* optional stuff */

  //do_what_you_want();

  return 1;
}


/* here the particle data is at your disposal 
 */
static void do_what_you_want(void){
    for(int i=1; i<=NumPart; i++)

    return;
}



/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
static void Gadget_unit_conversion(void){

    /* physical constants in cgs units */
    double GRAVITY   = GRAVITY_CONSTANT_CGS;
    double BOLTZMANN = BOLTZMANN_CONSTANT_CGS;
    double PROTONMASS = PROTON_MASS_CGS;

    /* internal unit system of the code */
    double UnitLength_in_cm= KPC_CGS;   /*  code length unit in cm/h */
    double UnitMass_in_g= 1e10*MSUN_CGS;         /*  code mass unit in g/h */
    double UnitVelocity_in_cm_per_s= VELOCITY_KMS_CGS;

    double UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    //double UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
    //double UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
    double UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

    //double G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);

    double Xh= 0.76;  /* mass fraction of hydrogen */
    double HubbleParam= 0.65;

    for(int i=1; i<=NumPart; i++){
        if(GadgetP[i].Type==0){  /* gas particle */
	        double MeanWeight= 4.0/(3*Xh+1+4*Xh*GadgetP[i].Ne) * PROTONMASS;

	        double u  = GadgetP[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;
	        double gamma= 5.0/3;
	        /* get temperature in Kelvin */
	        GadgetP[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	    }
    }

    return;
}



/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
static void Gadget_load_snapshot(char *fname, int files){
    FILE *fd;
    char   buf[MaxCharactersInLine];
    int    dummy,ntot_withmasses;
    int    pc,pc_new,pc_sph;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

    for(int i=0, pc=1; i<files; i++, pc=pc_new){
        if(files>1)
            sprintf(buf,"%s.%d",fname,i);
        else
            sprintf(buf,"%s",fname);

        FileOpen(fd,buf,"r");
        printf("reading `%s' ...\n",buf); fflush(stdout);

        fread(&dummy, sizeof(dummy), 1, fd);
        fread(&GadgetHeader1, sizeof(GadgetHeader1), 1, fd);
        fread(&dummy, sizeof(dummy), 1, fd);

        if(files==1){
            for(int k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                NumPart+= GadgetHeader1.npart[k];
            Ngas= GadgetHeader1.npart[0];
        }else{
            for(int k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                NumPart+= GadgetHeader1.npartTotal[k];
            Ngas= GadgetHeader1.npartTotal[0];
        }

        for(int k=0, ntot_withmasses=0; k<5; k++){
            if(GadgetHeader1.mass[k]==0)
                ntot_withmasses+= GadgetHeader1.npart[k];
        }

        if(i==0)
            Gadget_allocate_memory();

        SKIP;
        for(int k=0,pc_new=pc;k<6;k++){
            for(int n=0;n<GadgetHeader1.npart[k];n++){
                fread(&GadgetP[pc_new].Pos[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;

        SKIP;
        for(int k=0,pc_new=pc;k<6;k++){
            for(int n=0;n<GadgetHeader1.npart[k];n++){
                fread(&GadgetP[pc_new].Vel[0], sizeof(float), 3, fd);
                pc_new++;
            }
        }
        SKIP;
    

        SKIP;
        for(int k=0,pc_new=pc;k<6;k++){
            for(int n=0;n<GadgetHeader1.npart[k];n++){
                fread(&Id[pc_new], sizeof(int), 1, fd);
                pc_new++;
            }
        }
        SKIP;


        if(ntot_withmasses>0)
            SKIP;
        for(int k=0, pc_new=pc; k<6; k++){
            for(int n=0;n<GadgetHeader1.npart[k];n++){
                GadgetP[pc_new].Type=k;

                if(GadgetHeader1.mass[k]==0)
                    fread(&GadgetP[pc_new].Mass, sizeof(float), 1, fd);
                else
                    GadgetP[pc_new].Mass= GadgetHeader1.mass[k];

                pc_new++;
            }
        }
        if(ntot_withmasses>0)
            SKIP;
      

        if(GadgetHeader1.npart[0]>0){
            SKIP;
            for(int n=0, pc_sph=pc; n<GadgetHeader1.npart[0];n++){
                fread(&GadgetP[pc_sph].U, sizeof(float), 1, fd);
                pc_sph++;
            }
            SKIP;

            SKIP;
            for(int n=0, pc_sph=pc; n<GadgetHeader1.npart[0];n++){
                fread(&GadgetP[pc_sph].Rho, sizeof(float), 1, fd);
                pc_sph++;
            }
            SKIP;

            if(GadgetHeader1.flag_cooling){
                SKIP;
                for(int n=0, pc_sph=pc; n<GadgetHeader1.npart[0];n++){
                    fread(&GadgetP[pc_sph].Ne, sizeof(float), 1, fd);
                    pc_sph++;
                }
                SKIP;
            }else{
                for(int n=0, pc_sph=pc; n<GadgetHeader1.npart[0];n++){
                    GadgetP[pc_sph].Ne= 1.0;
                    pc_sph++;
                }
            }
        }
        fclose(fd);
    }


  Time= GadgetHeader1.time;
  Redshift= GadgetHeader1.time;
}


/* this routine allocates the memory for the 
 * particle data.
 */
static void Gadget_allocate_memory(void){

    fprintf(stderr,"allocating memory...\n");

    if(!(GadgetP=malloc(NumPart*sizeof(struct Gadget_particle_data)))){
        fprintf(stderr,"failed to allocate memory.\n");
        exit(0);
    }
    GadgetP--;   /* start with offset 1 */

    if(!(Id=malloc(NumPart*sizeof(int)))){
        fprintf(stderr,"failed to allocate memory.\n");
        exit(0);
    }
    Id--;   /* start with offset 1 */

    fprintf(stderr,"allocating memory...done\n");

    return;
}

/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
static void Gadget_reordering(void){

    fprintf(stderr,"reordering....\n");

    for(int i=1; i<=NumPart; i++){
        if(Id[i] != i){
            struct  Gadget_particle_data   psource= GadgetP[i];
            int idsource=Id[i];
            int dest=Id[i];

            do{
                struct  Gadget_particle_data psave= GadgetP[dest];
                int idsave=Id[dest];

                GadgetP[dest]= psource;
                Id[dest]= idsource;

                if(dest == i) 
                    break;

                psource= psave;
                idsource=idsave;

                dest=idsource;
            }while(1);
        }
    }

    fprintf(stderr,"done.\n");

    Id++;   
    free(Id);

    fprintf(stderr,"space for particle ID freed\n");
}

