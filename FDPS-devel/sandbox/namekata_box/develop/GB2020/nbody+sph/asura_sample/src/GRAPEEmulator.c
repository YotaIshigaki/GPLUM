#include "config.h"

#define EmuPipeLines    48

static int ni;
static int nj;
static double posi[EmuPipeLines][3];
static double epsi[EmuPipeLines];
static double kerneli[EmuPipeLines];
static double acceli[EmuPipeLines][3];
static double poti[EmuPipeLines];

static double (*posj)[3];
static double *mj;

static int Ninit = 0;
void InitializeGRAPEEmulator(void){
    
#if ( defined(HAVE_GRAPE7))||(defined(HAVE_GRAPE6A)||defined(HAVE_GRAPE5) )
    return;
#endif

#if !defined(JMEMSIZE)
#define JMEMSIZE 1  // dummy
//#error "JMEMSIZE is not defined."
//#error "Define JMEMSIZE to use the emulation mode."
#endif

    if(Ninit == 0){
        fprintf(stderr," ========================\n");
        fprintf(stderr," =GRAPE-5 emulation mode=\n");
        fprintf(stderr," ========================\n");
        posj = malloc(sizeof(double)*3*JMEMSIZE);
        mj = malloc(sizeof(double)*JMEMSIZE);
        Ninit ++;
    }


    return;
}

void g5_open_emu(void){
    return; // do nothing.
}

void g5_close_emu(void){
    return; // do nothing.
}

int g5_get_number_of_pipelines_emu(void){
    return EmuPipeLines;
}

int g5_get_number_of_boards_emu(void){
    return 1; // This function always returns 1.
}

void g5_set_range_emu(const double size1, const double size2, const double minimum_mass){
    return; // do nothing.
}

void g5_set_n_emu(const int number){
    nj = number;
    return;
}

void g5_set_xmj_emu(const int offset, const int number, 
        double Pos[restrict][3], double Mass[restrict]){
    nj = number;
    for(int i=0;i<number;i++){
        posj[i][0] = Pos[offset+i][0];
        posj[i][1] = Pos[offset+i][1];
        posj[i][2] = Pos[offset+i][2];
        mj[i] = Mass[offset+i];
    }

    return;
}

void g5_set_ip_emu(const int number_of_pipelines, double Pos[restrict][3],
        double Eps[restrict], double Hsml[restrict]){

    ni = number_of_pipelines;
    for(int i=0;i<number_of_pipelines;i++){
        posi[i][0] = Pos[i][0];
        posi[i][1] = Pos[i][1];
        posi[i][2] = Pos[i][2];
        epsi[i] = Eps[i];
        kerneli[i] = Hsml[i];
    }

    return;
}

void g5_run_emu(void){

    for(int i=0;i<ni;i++){
        acceli[i][0] = acceli[i][1] = acceli[i][2] = poti[i] = 0.e0;
        double eps2 = SQ(epsi[i]);
        for(int j=0;j<nj;j++){
            double dx[3] = {posi[i][0]-posj[j][0],
                            posi[i][1]-posj[j][1],
                            posi[i][2]-posj[j][2]}; 
            double r2 = NORM2(dx) + eps2;
            double r = sqrt(r2);
            double ir = 1.e0/r;
            double ir3 = CUBE(ir);
            acceli[i][0] -= mj[j]*ir3*dx[0];
            acceli[i][1] -= mj[j]*ir3*dx[1];
            acceli[i][2] -= mj[j]*ir3*dx[2];
            poti[i] -= mj[j]*ir;
        }
    }
    return;
}

void g5_get_force_emu(const int number_of_pipelines, 
        double Acc[restrict][3], double Pot[restrict]){

    for(int i=0;i<number_of_pipelines;i++){
        Acc[i][0] = acceli[i][0];
        Acc[i][1] = acceli[i][1];
        Acc[i][2] = acceli[i][2];
        Pot[i] = -poti[i];
    }
    return;
}

