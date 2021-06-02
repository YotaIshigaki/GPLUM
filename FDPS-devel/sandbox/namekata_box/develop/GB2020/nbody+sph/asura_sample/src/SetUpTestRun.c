#include "config.h"
#include "Logs.h"
#include "HydroMisc.h"
#include "Integral.h"
#include "PlantHydroTree.h"
#include "Decomposition.h"
#include "NeighborSearch.h"
#include "ForceMisc.h"
#include "ForceDirect.h"
#include "SetUpTestRun.h"
#include "HydroKernel.h"
#include "ThermalConductivity.h"
#include "StellarFeedback.h"
#include "EstimateXFe.h"


static void InitializeAllActiveParticleNumbers(void){

    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    long int Actives[5];
    Actives[0] = Pall.NActivesHydro;
    Actives[1] = Pall.NActivesStars;
    Actives[2] = Pall.NActivesSink;
    Actives[3] = Pall.NActivesDM;
    Actives[4] = Pall.NActives;
    long int GlobalActives[5];
    MPI_Allreduce(Actives,GlobalActives,5,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    Pall.NActivesHydro = GlobalActives[0];
    Pall.NActivesStars = GlobalActives[1];
    Pall.NActivesSink  = GlobalActives[2];
    Pall.NActivesDM    = GlobalActives[3];
    Pall.NActives      = GlobalActives[4];

    return ;
}

int FirstSetViscosityParameters = true;
static void SetConstantViscosityParameters(const double Alpha){

    // Pall.HydroAlpha = 1.0;
    Pall.HydroAlpha = Alpha;
    Pall.HydroBeta =  2.0*Pall.HydroAlpha;
    Pall.HydroEta2 =  0.01;

    FirstSetViscosityParameters = false;

    return;
}

static void SetVariableViscosityParameters(const double AlphaMin, const double AlphaMax, const double SourceFactor, const double InvDumpingScale){

    Pall.ViscousAlphaMin = AlphaMin; // 0.1
    Pall.ViscousAlphaMax = AlphaMax; // 
    Pall.ViscousS = SourceFactor;    // 1.0
    Pall.ViscousL = InvDumpingScale; // 0.1

    FirstSetViscosityParameters = false;

    return;
}

void SetViscosityParameters(const double AlphaMin, const double AlphaMax, const double SourceFactor, const double InvDumpingScale){
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(AlphaMin,AlphaMax,SourceFactor,InvDumpingScale);
#else // USE_VARIABLE_ALPHA //}//{
    double Alpha = AlphaMin;
    SetConstantViscosityParameters(Alpha);
#endif
    return;
}


void ActivateAllparticles(void){

    if(FirstSetViscosityParameters == true){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"Need to call SetVariableViscosityParameters before calling ActivateAllparticles.\n");
            fflush(NULL);
        }
        assert(FirstSetViscosityParameters == false);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
#ifdef USE_VARIABLE_ALPHA //{
        //Phydro[i]->Alpha = Pall.ViscousAlphaMax;
        Phydro[i]->Alpha = 1.0;
#endif // USE_VARIABLE_ALPHA //}
        Phydro[i]->Active = PhydroBody(i)->Active;
#ifdef USE_SPSPH //{
        // m/rho = Zw/y
        // y = Zw*rho/m
        // Z = m*y/rho
        // if(Phydro[i]->Rho > 0.e0){
            // Phydro[i]->PseudoDensity = Phydro[i]->Rho;
            // Phydro[i]->Zw = Phydro[i]->Mass;
        // } else {
            //Phydro[i]->Zw = 1.0;
            //Phydro[i]->PseudoDensity = Phydro[i]->Zw*Phydro[i]->Rho/Phydro[i]->Mass;
        // }
        if(Phydro[i]->Rho > 0.e0){
            Phydro[i]->PseudoDensity = 1.e0;
            Phydro[i]->Zw = Phydro[i]->PseudoDensity*Phydro[i]->Mass/Phydro[i]->Rho;
        } else {
            Phydro[i]->Zw = 1.0;
            Phydro[i]->PseudoDensity = Phydro[i]->Zw*Phydro[i]->Rho/Phydro[i]->Mass;
        }
#ifdef TASK_HYDROSTATIC //{
        // Phydro[i]->Zw = 1.0;
        // Phydro[i]->PseudoDensity = Phydro[i]->Zw*Phydro[i]->Rho/Phydro[i]->Mass;
#endif // TASK_HYDROSTATIC //}

#ifdef TASK_KELVINHELMHOLTZ_INSTABILITY //{
        Phydro[i]->PseudoDensity = 1.e0;
        Phydro[i]->Zw = Phydro[i]->PseudoDensity*Phydro[i]->Mass/Phydro[i]->Rho;
#endif // TASK_HYDROSTATIC //}

#endif //USE_SPSPH //}
    }

    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Active = ON;
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;
    }

    // count all active particles.
    int NActives = 0;
    int NActivesDM = 0;
    int NActivesHydro = 0;
    int NActivesStars = 0;
    int NActivesSink = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
            NActivesHydro ++;
        } else if(Pbody[i]->Type == TypeStar){
            NActivesStars ++;
        } else if(Pbody[i]->Type == TypeSink){
            NActivesSink ++;
        } else if(Pbody[i]->Type == TypeDM){
            NActivesDM ++;
        }
        NActives ++;
    }

    Pall.NActivesAll = Pall.NActives = NActives;
    Pall.NActivesDM = NActivesDM;
    Pall.NActivesHydro = NActivesHydro;
    Pall.NActivesStars = NActivesStars;
    Pall.NActivesSink = NActivesSink;

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    return ;
}



// 3D collapse problem.
void Init3DCollapseTest(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

    InitializeRandomGenerator(1977+MPIGetMyID());

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

	//double hmean = 2.e0/cbrt((double)number);

    // count Pall.Ntotal,Pall.Nhydro

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
                    
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}


    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = 1.e0/(double)count;
    double eps = 0.1*(cbrt((double)1736/(double)count));
    //double eps = 0.01*(cbrt((double)1736/(double)count));
    //eps = 0.01;
    double Uinit = 0.05;

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
				
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Pbody[mycount]->Vel[0] = TINY;
                        Pbody[mycount]->Vel[1] = TINY;
                        Pbody[mycount]->Vel[2] = TINY;

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 0.1*dx;
                        PbodyHydro(mycount)->U = Uinit;

                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    ActivateAllparticles();

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"3Dcollapse.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_DISPH //{
        double Smoothed = Phydro[i]->EnergyDensity;
        double Weight = Phydro[i]->Mass*Phydro[i]->U;
#elif defined(USE_SPSPH) //}//{
        double Smoothed = Phydro[i]->PseudoDensity;
        double Weight = Phydro[i]->Zw;
#else // USE_SPSPH //}//{
        double Smoothed = Phydro[i]->Rho;
        double Weight = Phydro[i]->Mass;
#endif // USE_DISPH //}
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i),
                Smoothed,Weight);
    }
    fclose(fp);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    //Pall.Ns = 32;
    Pall.Ns = 64;
    Pall.Npm = 2;

    Pall.TEnd = 3.0;
    //Pall.TEnd = 10.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/100.0; 
    //Pall.OutPutInterval = Pall.TEnd/1.0; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/3D.ASCII");
    strcpy(Pall.BaseFileName,"./data/3D");
    strcpy(Pall.RestartFileName,"./data/3D.dump");

    return;
}

#ifdef TASK_CLOUD_EXPLOSION
void InitCloudExplosion(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
                    
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double mass = 1.e0/(double)count;
    double eps = 0.1*(cbrt((double)1736/(double)count));
    double Uinit = 0.05;

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
				
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Pbody[mycount]->Vel[0] = TINY;
                        Pbody[mycount]->Vel[1] = TINY;
                        Pbody[mycount]->Vel[2] = TINY;

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 0.1*dx;
                        PbodyHydro(mycount)->U = Uinit;

                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    ReConnectPointers();
    UpdateTotalNumber();

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.TEnd = 10.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/1.0; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/EXP.ASCII");
    strcpy(Pall.BaseFileName,"./data/EXP");
    strcpy(Pall.RestartFileName,"./data/EXP.dump");
    InitializeRandomGenerator(1977+MPIGetMyID());

    return;
}
struct StructWriteCE{
    float Pos[3];
    float Vel[3];
    float Mass,Kernel,U;
};

void ReadInitCloudExplosion(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int Number;
    FILE *fp;
    FileOpen(fp,"./CE.init","rb");
    fread(&Number,sizeof(int),1,fp);
    dprintlmpi(Number);
    
    int mycount = 0;
	for(int i=0;i<Number;i++){
        if(i%NProcs == MyID){
            mycount ++;
        }
	}
    fclose(fp);
    dprintlmpi(Number);

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = Number;
    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double eps = 0.1*(cbrt((double)1736/(double)Number));
    for(int k=0;k<NProcs;k++){
        int mycount = 0;
        FileOpen(fp,"./CE.init","rb");
        int number;
        fread(&number,sizeof(int),1,fp);
        struct StructWriteCE WriteCE;
        dprintlmpi(number);
        for(int i=0;i<Number;i++){
            int ID;
            double Pos[3];
            double Kernel,Mass,U;

            fread(&WriteCE,sizeof(struct StructWriteCE),1,fp);
            if(i%NProcs == MyID){
                Pbody[mycount]->Active = ON;
                Pbody[mycount]->Use = ON;
                Pbody[mycount]->Type = TypeHydro;
                Pbody[mycount]->GlobalID = i; // read

                Pbody[mycount]->Pos[0] = WriteCE.Pos[0]; // read
                Pbody[mycount]->Pos[1] = WriteCE.Pos[1]; // read
                Pbody[mycount]->Pos[2] = WriteCE.Pos[2]; // read

                Pbody[mycount]->Vel[0] = WriteCE.Vel[0];
                Pbody[mycount]->Vel[1] = WriteCE.Vel[1];
                Pbody[mycount]->Vel[2] = WriteCE.Vel[2];

                Pbody[mycount]->Mass = WriteCE.Mass;        // read
                Pbody[mycount]->Eps = eps;

                PbodyHydro(mycount)->Use = ON;
                PbodyHydro(mycount)->Kernel = WriteCE.Kernel; // read
                PbodyHydro(mycount)->U = WriteCE.U;  // read
                mycount ++;
            }
        }
        fclose(fp);
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    ReConnectPointers();
    UpdateTotalNumber();

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.TEnd = 5.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = 1.0; 
    Pall.OutPutInterval = 0.2; 
    //Pall.OutPutInterval = 5.0; 
    //Pall.OutPutInterval = 0.1; 
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/EXP.ASCII");
    strcpy(Pall.BaseFileName,"./data/EXP");
    strcpy(Pall.RestartFileName,"./data/EXP.dump");
    InitializeRandomGenerator(1977+MPIGetMyID());

    return;
}

#endif

#ifdef TASK_BLAST_WAVE //{
void InitBlastWave(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = number;

    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
                if(count%NProcs == MyID){
                    mycount ++;
                }
                count ++;
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));

    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1.0;

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;
    Pall.BoxCenter[0] = Pall.BoxCenter[1] = Pall.BoxCenter[2] = 0.e0;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


	double dx = 1.e0/((double)NX);
	double Lhalf = 0.5;
    double mass = 1.e0/(double)count;
    double Uinit = 1.e-10/mass;
    int PeakIndex = 0;

	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
                double Pos[3];

				Pos[0] = i*dx-Lhalf;
				Pos[1] = j*dx-Lhalf;
				Pos[2] = k*dx-Lhalf;

                if(count%NProcs == MyID){
                    Pbody[mycount]->Active = ON;
                    Pbody[mycount]->Use = ON;
                    Pbody[mycount]->Type = TypeHydro;
                    Pbody[mycount]->GlobalID = count;

                    Pbody[mycount]->Pos[0] = Pos[0];
                    Pbody[mycount]->Pos[1] = Pos[1];
                    Pbody[mycount]->Pos[2] = Pos[2];

                    Pbody[mycount]->Vel[0] = TINY;
                    Pbody[mycount]->Vel[1] = 0.e0;
                    Pbody[mycount]->Vel[2] = 0.e0;

                    Pbody[mycount]->Mass = mass;
                    Pbody[mycount]->Eps = 1.0;

                    PbodyHydro(mycount)->Use = ON;
                    PbodyHydro(mycount)->Kernel = 2.0*dx;
                    PbodyHydro(mycount)->U = Uinit;
                    if( (i==NX/2)&&(j==NY/2)&&(k==NZ/2) ){
                        PeakIndex = mycount;
                        //PbodyHydro(mycount)->U = 1.e0/mass;
                    }

                    mycount ++;
                }
                count ++;
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }


    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    //Pall.Ns = 64;
    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = 0.041;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(3.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/BlastWave.ASCII");
    strcpy(Pall.BaseFileName,"./data/BlastWave");
    strcpy(Pall.RestartFileName,"./data/BlastWave.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 1.01;

    InitializeRandomGenerator(1977);


#if 0
    Phydro[PeakIndex]->U = 1.0/PhydroMass(PeakIndex);
#else
    //AddSmoothedEnergyBlastWave(PeakIndex);
    AddSmoothedEnergyBlastWaveDirectSearch(PeakIndex, number);
#endif


    return;
}

void ReadBlastWave(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    // file open and count.
    FILE *fp;
    char fname[MaxCharactersInLine];
    FileOpen(fp,"BlastWaveGlass.dat","r");
    int Ntotal = 0;
    while(fscanf(fp,"%*d %*g %*g %*g %*g %*g %*g") != EOF){
        Ntotal ++;
    }
    fclose(fp);

    int mycount = 0;
	for(int i=0;i<Ntotal;i++){
        if(i%NProcs == MyID){
            mycount ++;
        }
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,2.0,1.0,0.1);
#else
    SetConstantViscosityParameters(2.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.UnitTime = 1.0;
    Pall.UnitLength = 1.0;
    Pall.UnitMass = 1.0;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;
    Pall.BoxCenter[0] = Pall.BoxCenter[1] = Pall.BoxCenter[2] = 0.0;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = Ntotal;
    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    for(int k=0;k<NProcs;k++){
        int mycount = 0;
        FILE *fp;
        char fname[MaxCharactersInLine];
        sprintf(fname,"BlastWaveGlass.dat");
        FileOpen(fp,fname,"r");
        for(int i=0;i<Ntotal;i++){
            int ID;
            double Pos[3];
            double Kernel,Mass,U;

            // read line
            fscanf(fp,"%d %le %le %le %le %le %le",&ID,Pos,Pos+1,Pos+2,&Kernel,&U,&Mass);
            if(i%NProcs == MyID){
                //fprintf(stderr,"%d %le %le %le %le %le %le\n",ID,Pos[0],Pos[1],Pos[2],Kernel,U,Mass);
                Pbody[mycount]->Active = ON;
                Pbody[mycount]->Use = ON;
                Pbody[mycount]->Type = TypeHydro;
                //Pbody[mycount]->GlobalID = ID; // read
                Pbody[mycount]->GlobalID = i; 

                if(Pos[0] > Pall.Lboxh[0])
                    Pos[0] -= Pall.Lbox[0]; 
                if(Pos[1] > Pall.Lboxh[1])
                    Pos[1] -= Pall.Lbox[1]; 
                if(Pos[2] > Pall.Lboxh[2])
                    Pos[2] -= Pall.Lbox[2]; 

                if(Pos[0] < -Pall.Lboxh[0])
                    Pos[0] += Pall.Lbox[0]; 
                if(Pos[1] < -Pall.Lboxh[1])
                    Pos[1] += Pall.Lbox[1]; 
                if(Pos[2] < -Pall.Lboxh[2])
                    Pos[2] += Pall.Lbox[2]; 

                Pbody[mycount]->Pos[0] = Pos[0]; // read
                Pbody[mycount]->Pos[1] = Pos[1]; // read
                Pbody[mycount]->Pos[2] = Pos[2]; // read

                Pbody[mycount]->Vel[0] = 0.e0;
                Pbody[mycount]->Vel[1] = 0.e0;
                Pbody[mycount]->Vel[2] = 0.e0;

                Pbody[mycount]->Mass = Mass;        // read
                Pbody[mycount]->Eps = 1.0;

                PbodyHydro(mycount)->Use = ON;
                PbodyHydro(mycount)->Kernel = Kernel; // read
                PbodyHydro(mycount)->Rho = 1.0; // read
                PbodyHydro(mycount)->U = U;  // read
                mycount ++;
            }
        }
        fclose(fp);
	}

    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    //Pall.Ns = 64;
    Pall.Ns = 32;
    Pall.Npm = 2;
    // Pall.Ns = 128;
    // Pall.Npm = 8;

    //Pall.TEnd = 0.04;
    //Pall.TEnd = 0.081;
    //Pall.TEnd = 0.041;
    Pall.TEnd = 0.05;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/BlastWave.ASCII");
    strcpy(Pall.BaseFileName,"./data/BlastWave");
    strcpy(Pall.RestartFileName,"./data/BlastWave.dump");
    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = 0.005;
    //Pall.OutPutInterval = 0.01;
    //Pall.OutPutInterval = 0.02;
    Pall.OutPutInterval = 1.01;
    //Pall.OutPutInterval = 0.04;

    InitializeRandomGenerator(1977);

#if 1
#if 0
    Phydro[PeakIndex]->U = 1.0/PhydroMass(PeakIndex);
#else
    //AddSmoothedEnergyBlastWave(PeakIndex);
    int PeakIndex = 0;
    double d_min = NORM(Phydro[0]->PosP); 
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->U = 1.e-10/PhydroBody(i)->Mass;
        //Phydro[i]->U = 1.e-3/PhydroBody(i)->Mass;
        //Phydro[i]->U = 0.03*1.e-7;
        //PhydroBody(i)->Mass;
        if(d_min > NORM(Phydro[i]->PosP)){
            PeakIndex = i;
            d_min = NORM(Phydro[i]->PosP);
        }
    }
    dprintlmpi(PeakIndex);
    dprintlmpi((int)(cbrt(Pall.Nhydro)));
    AddSmoothedEnergyBlastWaveDirectSearch(PeakIndex, (int)(cbrt(Pall.Nhydro)));
    // search peak.
    double Upeak = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++)
        Upeak = fmax(Upeak,Phydro[i]->U);
    double Uambient = 1.1e-6*Upeak;
    //double Uambient = 1.e-5*Upeak;
    for(int i=0;i<Pall.Nhydro;i++)
        Phydro[i]->U += Uambient;
    fprintf(stderr,"Upeak and Uambient = %g and %g\n",Upeak,Uambient);
#endif
#endif

    return;
}

void InitBlastWaveRandom(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = number;

    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
                if(count%NProcs == MyID){
                    mycount ++;
                }
                count ++;
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = 0.5*Pall.Lbox[0];
    Pall.Lboxh[1] = 0.5*Pall.Lbox[1];
    Pall.Lboxh[2] = 0.5*Pall.Lbox[2];
    Pall.BoxCenter[0] = 0.e0;
    Pall.BoxCenter[1] = 0.e0;
    Pall.BoxCenter[2] = 0.e0;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    PbodySize = mycount;
    PhydroSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


	double dx = 1.e0/((double)NX);
    double mass = 1.e0/(double)count;
    /*
	double dx = 2.e0/((double)NX);
	double Lhalf = 1.0;
    double mass = 8.e0/(double)count;
    */
    double Uinit = 1.e-10/mass;

	for(int i=0;i<mycount;i++){
        double Pos[3];

        Pos[0] = gsl_rng_uniform(RandomGenerator)-0.5;
        Pos[1] = gsl_rng_uniform(RandomGenerator)-0.5;
        Pos[2] = gsl_rng_uniform(RandomGenerator)-0.5;

        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = mycount;

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        Pbody[i]->Pos[2] = Pos[2];

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = mass;
        Pbody[i]->Eps = 1.0;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = 2.0*dx;
        PbodyHydro(i)->U = Uinit;
    }
    
#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"InitBlastWave.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    exit(0);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    //Pall.Npm = 0;

    Pall.TEnd = 0.081;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    strcpy(Pall.ASCIIFileName,"./data/BlastWave.ASCII");
    strcpy(Pall.BaseFileName,"./data/BlastWave");
    strcpy(Pall.RestartFileName,"./data/BlastWave.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.01;

    InitLogFiles();

    return;
}

static double BlastWaveKernel(const double r, const double InvKerneli){

    double u = r*InvKerneli;
    const static double coef3d = M_1_PI;
    //double coef = coef3d*CUBE(InvKerneli);
    double coef = 1.e0;

    if(u<1.e0){
        return (coef*(1.e0 - (1.5)*SQ(u) + (0.75)*CUBE(u)));
    } else if (u<2.e0){
        return (coef*((0.25)*CUBE(2.e0-(u))));
    } else {
        return (0.e0);
    }
}

struct StructDistance{
    int index;
    double dist;
};

static int DistCmp(const void *x, const void *y){
    const struct StructDistance *pointer1 = x;
    const struct StructDistance *pointer2 = y;
    if( pointer1->dist > pointer2->dist)
        return 1;
    else if( pointer1->dist < pointer2->dist)
        return -1;
    else
        return 0;
}

//#define BlastSmoothedNumber  (128)
#define BlastSmoothedNumber  (32)
//#define BlastSmoothedNumber  (4)
void AddSmoothedEnergyBlastWaveDirectSearch(const int PeakIndex, const int number){

    if(MPIGetNumProcs()>1)
        exit(-1);
    
    struct StructDistance *Distances; 
    Distances = malloc(sizeof(struct StructDistance)*Pall.Nhydro);
    for(int i=0;i<Pall.Nhydro;i++){
        //Distances[i].dist = NORM(Phydro[i]->PosP);
        Distances[i].dist = DISTANCE(Phydro[i]->PosP,Phydro[PeakIndex]->PosP);
        Distances[i].index = i;
    }

    qsort(Distances,Pall.Nhydro,sizeof(struct StructDistance),(int(*)(const void*, const void*))DistCmp);

    // Peak first Pall.Ns particles.
    int Neighbors[BlastSmoothedNumber];
    for(int i=0;i<BlastSmoothedNumber;i++){
        Neighbors[i] = Distances[i].index;
    }
    double kernel = 0.5*DISTANCE(Phydro[Neighbors[0]]->PosP,Phydro[Neighbors[BlastSmoothedNumber-1]]->PosP);

    double wt = 0;
    for(int i=0;i<BlastSmoothedNumber;i++){
        double r = DISTANCE(Phydro[PeakIndex]->PosP,Phydro[Neighbors[i]]->PosP);
        //double w = BlastWaveKernel(r,1.e0/Phydro[PeakIndex]->Kernel);
        double w = BlastWaveKernel(r,1.e0/kernel);
        wt += w;
    }

    FILE *fp;
    FileOpen(fp,"heated.dat","w");
    fprintf(fp,"%d\n",BlastSmoothedNumber);

    double iwt = 1.e0/wt;
    eprintlmpi(iwt);
    double Uinit = 1.0/Phydro[PeakIndex]->Mass;
    double totale = 0.0;
    for(int i=0;i<BlastSmoothedNumber;i++){
        double r = DISTANCE(Phydro[PeakIndex]->PosP,Phydro[Neighbors[i]]->PosP);
        //double w = BlastWaveKernel(r,1.e0/Phydro[PeakIndex]->Kernel);
        double w = BlastWaveKernel(r,1.e0/kernel);
        fprintf(stderr,"Add energy[%d] %g -> %g, dU = %g, r %g r/h %g \n",i,
                //Phydro[Neighbors[i]]->U,Phydro[Neighbors[i]]->U+Uinit*w*iwt,Uinit*w*iwt,r,r/Phydro[PeakIndex]->Kernel);
                Phydro[Neighbors[i]]->U,Phydro[Neighbors[i]]->U+Uinit*w*iwt,Uinit*w*iwt,r,r/kernel);
        Phydro[Neighbors[i]]->U += Uinit*w*iwt;
        fprintf(fp,"%ld %g %g %g\n",PhydroBody(Neighbors[i])->GlobalID,
                Phydro[Neighbors[i]]->PosP[0],
                Phydro[Neighbors[i]]->PosP[1],
                Phydro[Neighbors[i]]->PosP[2]);
        totale += Uinit*w*iwt;
    }
    gprintlmpi(totale);
    fclose(fp);
    fflush(fp);

#if 1 //cheker
    FILE *fh;
    char fname[MaxCharactersInLine];
    sprintf(fname,"heated.dat");
    FileOpen(fh,fname,"r");
    int Nsmooth;
    fscanf(fh,"%d",&Nsmooth);
    dprintlmpi(Nsmooth);
    int GID[Nsmooth];
    for(int i=0;i<Nsmooth;i++){
        fscanf(fh,"%d %*g %*g %*g",GID+i);
        fprintf(stderr,"GDI = %d\n",GID[i]);
    }
    fclose(fh);

    sprintf(fname,"heated_check.dat");
    FileOpen(fh,fname,"w");
    for(int k=0;k<Nsmooth;k++){
        for(int i=0;i<Pall.Nhydro;i++){
            if(GID[k] == PhydroBody(i)->GlobalID){
                fprintf(fh,"%ld %g %g %g\n",PhydroBody(i)->GlobalID,
                        Phydro[i]->PosP[0],
                        Phydro[i]->PosP[1],
                        Phydro[i]->PosP[2]);
            }
        }
    }
    fclose(fh);
#endif
    free(Distances);

    return;
}

void AddSmoothedEnergyBlastWave(const int number){

    int NX,NY,NZ;
    NX = NY = NZ = number;

    if(MPIGetNumProcs()>1)
        exit(-1);

    int count = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
                if((i==NX/2)&&(j==NY/2)&&(k==NZ/2)){
                    int Neighbors[MaxNeighborSize];
                    int Nlist = GetNeighborsLimited(Pbody[count]->Pos,2.0*PbodyHydro(count)->Kernel,Neighbors);
                    double wt = 0;
                    for(int leaf=0;leaf<Nlist;leaf++){
                        double r = DISTANCE(Pbody[count]->Pos,Pbody[Neighbors[leaf]]->Pos);
                        double w = BlastWaveKernel(r,1.e0/PbodyHydroKernel(count));
                        wt += w;
                    }
                    double iwt = 1.e0/wt;
                    eprintlmpi(iwt);
                    double Uinit = 1.0/Pbody[count]->Mass;
                    double totale = 0; 
                    for(int leaf=0;leaf<Nlist;leaf++){
                        double r = DISTANCE(Pbody[count]->Pos,Pbody[Neighbors[leaf]]->Pos);
                        double w = BlastWaveKernel(r,1.e0/PbodyHydroKernel(count));
                        //PbodyHydro(Neighbors[leaf])->U = w*iwt;
                        //PbodyHydro(Neighbors[leaf])->U = Uinit*w*iwt;
                        if(count == leaf){
                            PbodyHydro(Neighbors[leaf])->U = Uinit*w*iwt;
                            totale += Uinit*w*iwt;
                        }else if(PbodyHydro(Neighbors[leaf])->U < Uinit*w*iwt){
                            PbodyHydro(Neighbors[leaf])->U = Uinit*w*iwt;
                            totale += Uinit*w*iwt;
                        }
                    }
                    eprintlmpi(totale);
                    eprintlmpi(totale*Pbody[count]->Mass);
                }
                count ++;
			}
		}
	}
    return;
}

void GlassCondition(const double RhoTrue){

    if(MPIGetNumProcs()>1)
        exit(-1);

    double iRhoTrue = 1.e0/RhoTrue;
    double epsilon = 1.e0;
    double kernel = cbrt((3.0*Pall.Ns/(32.0*M_PI*Pall.Nhydro))*Pall.Lbox[0]*Pall.Lbox[1]*Pall.Lbox[2]);

    double *RhoOld;
    RhoOld = malloc(sizeof(double)*Pall.Nhydro);

    double (*PosNew)[3];
    PosNew = malloc(sizeof(double)*3*Pall.Nhydro);
    
    int iteration = 0;
    do{
        // get mean and max;
        double RhoMean = 0.e0;
        for(int i=0;i<Pall.Nhydro;i++){
            RhoMean += PbodyHydroRho(i);
        }
        RhoMean /= Pall.Nhydro;

        double RhoMax = 0.e0;
        for(int i=0;i<Pall.Nhydro;i++){
            RhoMax = fmax(PbodyHydroRho(i),RhoMax);
        }

        double RhoVar = 0.e0;
        for(int i=0;i<Pall.Nhydro;i++){
            RhoVar += SQ(PbodyHydroRho(i)-RhoMean);
        }
        RhoVar/=Pall.Nhydro;
        RhoVar = sqrt(RhoVar);

        fprintf(stderr,"Iteration [%d]",iteration);
        fprintf(stderr," :Mean density = %g, the variance = %g, Max Density = %g\n",RhoMean,RhoVar,RhoMax);
        if((1.01*RhoTrue>(RhoMean+RhoVar))&&(0.99*RhoTrue<(RhoMean+RhoVar)))
            break;

        for(int i=0;i<Pall.Nhydro;i++){
            RhoOld[i] = PbodyHydroRho(i);
            PosNew[i][0] = PosNew[i][1] = PosNew[i][2] = 0.e0;  
        }

        int Neighbors[MaxNeighborSize];
        for(int i=0;i<Pall.Nhydro;i++){
            int Nlist = GetNeighborsLimited(Pbody[i]->Pos,2.0*kernel,Neighbors);
            PbodyHydroRho(i) = 0.e0;
            for(int k=0;k<Nlist;k++){
                int leaf = Neighbors[k];
                double Posj[3] = {Pbody[leaf]->Pos[0],Pbody[leaf]->Pos[1],Pbody[leaf]->Pos[2]};
                double xij[3];
#ifdef PERIODIC_RUN 
                for(int l=0;l<DIMENSION;l++){
                    xij[l] = PeriodicDistance(Pbody[i]->Pos[l],Posj[l],l);
                }
#else // PERIODIC_RUN
                for(int l=0;l<DIMENSION;l++){
                    xij[l] = Pbody[i]->Pos[l]-Posj[l];
                }
#endif // PERIODIC_RUN

                double r = NORM(xij);
                double w = 0.e0;
                double u = r/kernel;

                if(u<1.e0){
                    w = M_1_PI*(1.e0 - (1.5)*SQ(u) + (0.75)*CUBE(u));
                } else if (u<2.e0){
                    w = M_1_PI*((0.25)*CUBE(2.e0-(u)));
                } else {
                    w = 0.e0;
                }
                PbodyHydroRho(i) += Pbody[leaf]->Mass*w/CUBE(kernel);
                if(r > 0.e0){
                    double fact = epsilon*kernel*((RhoOld[leaf]-RhoTrue)*iRhoTrue)*w/r;
                    PosNew[i][0] += fact*xij[0];
                    PosNew[i][1] += fact*xij[1];
                    PosNew[i][2] += fact*xij[2];
                }
            }
        }

        for(int i=0;i<Pall.Nhydro;i++){
            Pbody[i]->Pos[0] += PosNew[i][0];
            Pbody[i]->Pos[1] += PosNew[i][1];
            Pbody[i]->Pos[2] += PosNew[i][2];
        }

        PeriodicWrapping();
        BuildPredictors();
        PlantHydroTree();

        FILE *fp;
        char fname[MaxCharactersInLine];
        sprintf(fname,"MakeGlass.%02d",MPIGetMyID());
        FileOpen(fp,fname,"w");
        for(int i=0;i<Pall.Nhydro;i++){
            fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                    PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                    Phydro[i]->Kernel,Phydro[i]->Rho,PhydroMass(i));
        }
        fclose(fp);

        iteration ++;
    }while(iteration<500);
    
    free(RhoOld);
    free(PosNew);

    double r2 = NORM(PhydroPos(0));
    int index = 0; 
    for(int i=0;i<Pall.Nhydro;i++){
        double d2 = NORM2(PhydroPos(i));
        if(r2 > d2){
            r2 = d2;
            index = i; 
        }
    }
    Phydro[index]->U = 1.e0/PhydroMass(index);

    return;
}
#endif // TASK_BLAST_WAVE //}

#ifdef TASK_SINUSOIDAL_WAVE //{
void InitSinusoidalWave(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.GravConst = 
    Pall.ConvertTtoU = Pall.ConvertUtoT = 
    Pall.ConvertDensityToCGS = Pall.ConvertNumberDensityToCGS = 1.0;

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    int NParticles = 100; 
    int AllocationSize = NParticles;
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double m = 1.0/(double)NParticles;
    double dx = 1.0/(double)NParticles;
    double rho = 1.0;
    double P = 1.0;
    double u = 1.0/Pall.Gm1;
    double cs = sqrt(Pall.GGm1*rho*u);
    double Vfact = 1.e-4;

    for(int i=0;i<NParticles;i++){
        PhydroBody(i)->Active = ON;
        PhydroBody(i)->Use = ON;
        PhydroBody(i)->Type = TypeHydro;
        PhydroBody(i)->GlobalID = i;

        double x = dx*(i+0.5);

        Phydro[i]->Kernel = 2*dx;
        PhydroBody(i)->Pos[0] = x;
        PhydroBody(i)->Pos[1] = PhydroBody(i)->Pos[2] = 0.0;

        double rad = 2*M_PI*x;
        //fprintf(stderr,"%g %g %g\n",x,rad,sin(rad));
        PhydroBody(i)->Vel[0] = Vfact*cs*sin(rad);
        PhydroBody(i)->Vel[1] = PhydroBody(i)->Vel[2] = 0.0;

        PhydroBody(i)->Mass = Phydro[i]->Mass = m;
        Phydro[i]->U = u;
        Phydro[i]->Rho = rho;

        Phydro[i]->Du =
        Phydro[i]->HydroAcc[0] =
        Phydro[i]->HydroAcc[1] =
        Phydro[i]->HydroAcc[2] = 0.e0;

        Phydro[i]->Use = ON; 
        PhydroBody(i)->Eps = 1.0;
#ifdef USE_PARTICLE_TAG
        Phydro[i]->Tag = 1; 
#endif // USE_PARTICLE_TAG
    }


    Pall.Ntotal = Pall.Nhydro = 
    Pall.Ntotal_t = Pall.Nhydro_t = NParticles;
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
        Pbody[i]->Active = ON;
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;
    }

    double TSound = 1.0/cs;
    double TEnd = 10*TSound;

    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Sound crossing time is %g, Tend = %g\n",TSound,TEnd);

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 5;
    Pall.Npm = 2;

    Pall.TEnd = TEnd;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = Pall.TEnd/50;
    Pall.OutPutInterval = Pall.TEnd/100;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/Sin.ASCII");
    strcpy(Pall.BaseFileName,"./data/Sin");
    strcpy(Pall.RestartFileName,"./data/Sin.dump");

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    return ;
}
#endif //TASK_SINUSOIDAL_WAVE //{

void InitSelfSimilarCooling(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    PbodySize = mycount;
    PhydroSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    // Set Unit;
    Pall.UnitLength  = 1.e0;
    Pall.UnitMass  = 1.e0;
    Pall.UnitTime = 1.e0;
    Pall.GravConst  = 2.0/(9*51.0);


    double mass = 1.e0/(double)count;
    double eps = 0.1;
    double Uinit = 0.05;
    double RhoH = 3.0/(4.0*M_PI);
    double Hi = sqrt(2.e0*Pall.GravConst);


	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Pbody[mycount]->Vel[0] = Hi*Pos[0];
                        Pbody[mycount]->Vel[1] = Hi*Pos[1];
                        Pbody[mycount]->Vel[2] = Hi*Pos[2];

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 2.0*dx;
                        PbodyHydro(mycount)->U = Uinit;

                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"SelfSimilarCooling.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = 61.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    InitLogFiles();

    return;
}

void ExternalPotentialForSelfSimilarCooling(void){

    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroActive(i)){
            double r2 = NORM2(PhydroPosP(i));
            double R2 = r2+SQ(0.1);
            double Ffact = Pall.GravConst*0.05/(R2*sqrt(R2));
            PhydroBody(i)->Acc[0] -= Ffact*PhydroPosP(i)[0];
            PhydroBody(i)->Acc[1] -= Ffact*PhydroPosP(i)[1];
            PhydroBody(i)->Acc[2] -= Ffact*PhydroPosP(i)[2];
        }
    }

    return; 
}

void Init3DShockTubeGrid(const int number){

    //int MyID = MPIGetMyID();
    //int NProcs = MPIGetNumProcs();

    if(MPIGetNumProcs()>1){
        fprintf(stderr,"This routine is not supported a parallel run.");
        exit(-1);
    }
    if((number%5 != 0)||((number/5)%16 != 0)){
        fprintf(stderr,"Particle number is not valid.");
        exit(-1);
    }

#if DIMENSION == 1
    fprintf(stderr,"The DIMENSION is not 3!");
    fprintf(stderr,"Check your config.h and have to set DIMENSION into 3!");
#endif

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));

    Pall.Lbox[0] = 1.e0;
    Pall.Lbox[1] = 1.e0/16.0;
    Pall.Lbox[2] = 1.e0/16.0;
    Pall.Lboxh[0] = 0.5*Pall.Lbox[0];
    Pall.Lboxh[1] = 0.5*Pall.Lbox[1];
    Pall.Lboxh[2] = 0.5*Pall.Lbox[2];
    Pall.BoxCenter[0] = 0.e0;
    Pall.BoxCenter[1] = 0.e0;
    Pall.BoxCenter[2] = 0.e0;

    Pall.Ntotal = Pall.Nhydro = number;
    Pall.Ntotal_t = Pall.Nhydro_t = number;
    int AllocationSize = number; 
    PbodySize = number;
    PhydroSize = number;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = SQ(1.0/8.0)/(double)number;

    int count = 0;

    int Nleft = (4*number/5);
    double dxleft = cbrt(1.0/(64*Nleft));
    int NleftBase = (1.0/8.0)/dxleft;
    // left
    // -1<x<0,-1/16<y<1/16,-1/16<z<1/16
    dprintlmpi(Nleft);
    dprintlmpi(NleftBase);
    eprintlmpi(dxleft);

	for(int i=0;i<NleftBase;i++){
	for(int j=0;j<NleftBase;j++){
	for(int k=0;k<8*NleftBase;k++){
        Pbody[count]->Pos[0] = dxleft*k-1.0;
        Pbody[count]->Pos[1] = dxleft*j-1.0/16.0;
        Pbody[count]->Pos[2] = dxleft*i-1.0/16.0;

        Pbody[count]->Active = ON;
        Pbody[count]->Use = ON;
        Pbody[count]->Type = TypeHydro;
        Pbody[count]->GlobalID = count;

        Pbody[count]->Vel[0] = TINY;
        Pbody[count]->Vel[1] = TINY;
        Pbody[count]->Vel[2] = TINY;

        Pbody[count]->Mass = mass;
        Pbody[count]->Eps = 1.0;

        PbodyHydro(count)->Use = ON;
        PbodyHydro(count)->Kernel = 2.0*dxleft;
        PbodyHydro(count)->U = 2.5;
        count ++;
    }
    }
    }


    int Nright = number-Nleft;
    double dxright = cbrt(1.0/(64*Nright));
    int NrightBase = (1.0/8.0)/dxright;
    // right
    // 0<x<1,-1/16<y<1/16,-1/16<z<1/16
    dprintlmpi(Nright);
    dprintlmpi(NrightBase);
    eprintlmpi(dxright);
	for(int i=0;i<NrightBase;i++){
	for(int j=0;j<NrightBase;j++){
	for(int k=0;k<8*NrightBase;k++){
        Pbody[count]->Pos[0] = dxright*k;
        Pbody[count]->Pos[1] = dxright*j-1.0/16.0;
        Pbody[count]->Pos[2] = dxright*i-1.0/16.0;

        Pbody[count]->Active = ON;
        Pbody[count]->Use = ON;
        Pbody[count]->Type = TypeHydro;
        Pbody[count]->GlobalID = count;

        Pbody[count]->Vel[0] = TINY;
        Pbody[count]->Vel[1] = TINY;
        Pbody[count]->Vel[2] = TINY;

        Pbody[count]->Mass = mass;
        Pbody[count]->Eps = 1.0;

        PbodyHydro(count)->Use = ON;
        PbodyHydro(count)->Kernel = 2.0*dxright;
        PbodyHydro(count)->U = 1.795;
        count ++;
    }
    }
    }
    dprintlmpi(count);
    dprintlmpi(number);
    assert(count == number);
    
#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"Init3DShockTube.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    exit(0);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    //Pall.Npm = 0;

    Pall.TEnd = 0.15;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    strcpy(Pall.ASCIIFileName,"./data/3DShockTube.ASCII");
    strcpy(Pall.BaseFileName,"./data/3DShockTube");
    strcpy(Pall.RestartFileName,"./data/3DShockTube.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1;

    InitLogFiles();
    InitializeRandomGenerator(1977);

    return;
}

void Init3DShockTube(const int number){

    //int MyID = MPIGetMyID();
    //int NProcs = MPIGetNumProcs();

    if(MPIGetNumProcs()>1){
        fprintf(stderr,"This routine is not supported a parallel run.");
        exit(-1);
    }
#if DIMENSION == 1
    fprintf(stderr,"The DIMENSION is not 3!");
    fprintf(stderr,"Check your config.h and have to set DIMENSION into 3!");
#endif

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.Lbox[0] = 1.e0;
    Pall.Lbox[1] = 1.e0/16.0;
    Pall.Lbox[2] = 1.e0/16.0;
    Pall.Lboxh[0] = 0.5*Pall.Lbox[0];
    Pall.Lboxh[1] = 0.5*Pall.Lbox[1];
    Pall.Lboxh[2] = 0.5*Pall.Lbox[2];
    Pall.BoxCenter[0] = 0.e0;
    Pall.BoxCenter[1] = 0.e0;
    Pall.BoxCenter[2] = 0.e0;

    Pall.Ntotal = Pall.Nhydro = number;
    Pall.Ntotal_t = Pall.Nhydro_t = number;
    int AllocationSize = number; 
    PbodySize = number;
    PhydroSize = number;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = SQ(1.0/8.0)/(double)number;
    double dv =  SQ(1.0/8.0);

    int Nleft = (4*number/5);
    double dxleft = cbrt(dv/Nleft);
    // left
    // -1<x<0,-1/16<y<1/16,-1/16<z<1/16
	for(int i=0;i<Nleft;i++){
        double Pos[3];

        Pos[0] = -1.0*gsl_rng_uniform(RandomGenerator);
        Pos[1] = (1.0/16.0)*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        Pos[2] = (1.0/16.0)*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);

        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = i;

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        Pbody[i]->Pos[2] = Pos[2];

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = mass;
        Pbody[i]->Eps = 1.0;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = 2.0*dxleft;
        PbodyHydro(i)->U = 2.5;
    }

    int Nright = number-Nleft;
    double dxright = cbrt(dv/Nright);
    // right
    // 0<x<1,-1/16<y<1/16,-1/16<z<1/16
	for(int i=Nleft;i<number;i++){
        double Pos[3];

        Pos[0] = 1.0*gsl_rng_uniform(RandomGenerator);
        Pos[1] = (1.0/16.0)*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        Pos[2] = (1.0/16.0)*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);

        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = i;

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        Pbody[i]->Pos[2] = Pos[2];

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = mass;
        Pbody[i]->Eps = 1.0;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = 2.0*dxright;
        PbodyHydro(i)->U = 1.795;
    }
    
#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"Init3DShockTube.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    //Pall.Npm = 0;

    Pall.TEnd = 0.15;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    strcpy(Pall.ASCIIFileName,"./data/3DShockTube.ASCII");
    strcpy(Pall.BaseFileName,"./data/3DShockTube");
    strcpy(Pall.RestartFileName,"./data/3DShockTube.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1;

    InitLogFiles();

    return;
}

void Read3DShockTube(const int number, char LoadFileName[]){

    int MyID = MPIGetMyID();

    if(MPIGetNumProcs()>1){
        fprintf(stderr,"This routine is not supported a parallel run.");
        exit(-1);
    }

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.Lbox[0] = 2.e0;
    Pall.Lbox[1] = 1.e0/16.0;
    Pall.Lbox[2] = 1.e0/16.0;
    Pall.Lboxh[0] = 0.5*Pall.Lbox[0];
    Pall.Lboxh[1] = 0.5*Pall.Lbox[1];
    Pall.Lboxh[2] = 0.5*Pall.Lbox[2];
    Pall.BoxCenter[0] = 0.e0;
    Pall.BoxCenter[1] = 0.e0;
    Pall.BoxCenter[2] = 0.e0;

    Pall.Ntotal = Pall.Nhydro = number;
    Pall.Ntotal_t = Pall.Nhydro_t = number;
    int AllocationSize = number; 
    PbodySize = number;
    PhydroSize = number;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    // read;
    FILE *fp;
    char fname[MaxCharactersInLine];
    FileOpen(fp,LoadFileName,"r");
	for(int i=0;i<number;i++){
        int ID;
        double Pos[3];
        double U,Kernel,Mass;
        
        fscanf(fp,"%d %le %le %le %*e %*e %*e %le %le %le %*d",
                &ID,Pos,Pos+1,Pos+2,&Kernel,&U,&Mass);
        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;

        Pbody[i]->GlobalID = ID; //read;

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        Pbody[i]->Pos[2] = Pos[2];

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = Mass; // read;
        Pbody[i]->Eps = 1.0;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = Kernel;  //read;
        PbodyHydro(i)->U = U; //read;

    }
    fclose(fp);

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

#if 1
    sprintf(fname,"Init3DShockTube.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = 0.15;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    strcpy(Pall.ASCIIFileName,"./data/3DShockTube.ASCII");
    strcpy(Pall.BaseFileName,"./data/3DShockTube");
    strcpy(Pall.RestartFileName,"./data/3DShockTube.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1;

    InitLogFiles();

    return;
}

void GlassCondition3DShockTube(void){

    if(MPIGetNumProcs()>1)
        exit(-1);


    double epsilon = 1.e0;
    double kernel = cbrt(3.0*Pall.Ns/(32.0*M_PI*Pall.Nhydro))*1.e0;

    double *RhoOld;
    RhoOld = malloc(sizeof(double)*Pall.Nhydro);

    double (*PosNew)[3];
    PosNew = malloc(sizeof(double)*3*Pall.Nhydro);
    
    int iteration = 0;
    do{
        // get mean and max;
        double RhoMean = 0.e0;
        for(int i=0;i<Pall.Nhydro;i++){
            RhoMean += PbodyHydroRho(i);
        }
        RhoMean /= Pall.Nhydro;

        double RhoMax = 0.e0;
        for(int i=0;i<Pall.Nhydro;i++){
            RhoMax = fmax(PbodyHydroRho(i),RhoMax);
        }

        double RhoVar = 0.e0;
        for(int i=0;i<Pall.Nhydro;i++){
            RhoVar += SQ(PbodyHydroRho(i)-RhoMean);
        }
        RhoVar/=Pall.Nhydro;
        RhoVar = sqrt(RhoVar);

        fprintf(stderr,"Iteration [%d]",iteration);
        fprintf(stderr," :Mean density = %g, the variance = %g, Max Density = %g\n",RhoMean,RhoVar,RhoMax);
        //if((1.01>(RhoMean+RhoVar))&&(0.99<(RhoMean+RhoVar)))
        if(fabs(RhoVar)<0.01)
            break;

        for(int i=0;i<Pall.Nhydro;i++){
            RhoOld[i] = PbodyHydroRho(i);
            PosNew[i][0] = PosNew[i][1] = PosNew[i][2] = 0.e0;  
        }

        int Neighbors[MaxNeighborSize];
        for(int i=0;i<Pall.Nhydro;i++){
            int Nlist = GetNeighborsLimited(Pbody[i]->Pos,2.0*kernel,Neighbors);
            PbodyHydroRho(i) = 0.e0;
            for(int k=0;k<Nlist;k++){
                int leaf = Neighbors[k];
                double Posj[3] = {Pbody[leaf]->Pos[0],Pbody[leaf]->Pos[1],Pbody[leaf]->Pos[2]};
                double xij[3];
#ifdef PERIODIC_RUN 
                for(int l=0;l<DIMENSION;l++){
                    xij[l] = PeriodicDistance(Pbody[i]->Pos[l],Posj[l],l);
                }
#else // PERIODIC_RUN
                for(int l=0;l<DIMENSION;l++){
                    xij[l] = Pbody[i]->Pos[l]-Posj[l];
                }
#endif // PERIODIC_RUN
                double r = NORM(xij);
                double w = 0.e0;
                double u = r/kernel;

                if(u<1.e0){
                    w = M_1_PI*(1.e0 - (1.5)*SQ(u) + (0.75)*CUBE(u));
                } else if (u<2.e0){
                    w = M_1_PI*((0.25)*CUBE(2.e0-(u)));
                } else {
                    w = 0.e0;
                }
                PbodyHydroRho(i) += Pbody[leaf]->Mass*w/CUBE(kernel);
                if(r > 0.e0){
                    double fact;
                    if(Pbody[i]->Pos[0] < 0.e0){
                        fact = epsilon*kernel*((RhoOld[leaf]-1.e0)/1.e0)*w/r;
                    } else {
                        //fact = epsilon*kernel*((RhoOld[leaf]-0.25)/0.25)*w/r;
                        fact = epsilon*kernel*((RhoOld[leaf]-1.e0)/1.e0)*w/r;
                    }
                    PosNew[i][0] += fact*xij[0];
                    PosNew[i][1] += fact*xij[1];
                    PosNew[i][2] += fact*xij[2];
                }
            }
        }

        for(int i=0;i<Pall.Nhydro;i++){
            Pbody[i]->Pos[0] += PosNew[i][0];
            Pbody[i]->Pos[1] += PosNew[i][1];
            Pbody[i]->Pos[2] += PosNew[i][2];
        }

        PeriodicWrapping();
        BuildPredictors();
        PlantHydroTree();

        FILE *fp;
        char fname[MaxCharactersInLine];
        sprintf(fname,"Init3DShockTube.%02d",MPIGetMyID());
        FileOpen(fp,fname,"w");
        for(int i=0;i<Pall.Nhydro;i++){
            fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                    PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                    Phydro[i]->Kernel,Phydro[i]->Rho,PhydroMass(i));
        }
        fclose(fp);

        iteration ++;
    }while(iteration<500);
    
    free(RhoOld);
    free(PosNew);


    return;
}

void OutPut3DCollapseTest(char *name){

	int n;
	double r,r2,dr,v_r,rho_r,u_r,p_r,cs_r;
    char fname[MaxCharactersInLine];
    FILE *fp;

    sprintf(fname,"%s.point.%03d.%03d",name,MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
		//r = NORM(PhydroPos(i));
        double Pos[3] = {PeriodicDistance(PhydroPos(i)[0],0,0),
                         PeriodicDistance(PhydroPos(i)[1],0,1),
                         PeriodicDistance(PhydroPos(i)[2],0,2)};
		r = NORM(Pos);
		if( r > 0 ){
			v_r = DOT_PRODUCT(PhydroPos(i),PhydroVel(i))/r;
#ifdef USE_DISPH //{
            double SmoothedValue = Phydro[i]->EnergyDensity;
            double Weight = Phydro[i]->Mass*Phydro[i]->U;
#elif defined(USE_SPSPH) //}//{
            double SmoothedValue = Phydro[i]->PseudoDensity;
            double Weight = Phydro[i]->Zw;
#else //}//{
            double SmoothedValue = Phydro[i]->Rho;
            double Weight = Phydro[i]->Mass;
#endif // USE_DISPH //}
#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
            double SmoothedNumber = Phydro[i]->SmoothedNumber;
#else
            double SmoothedNumber = Phydro[i]->Nlist;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}
#ifdef USE_VARIABLE_ALPHA //{
            double Alpha = Phydro[i]->Alpha;
            double DAlpha = Phydro[i]->DAlpha;
#else
            double Alpha = Pall.HydroAlpha;
            double DAlpha = 0.e0;
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //}

			fprintf(fp,"%ld %g %g %g %g %g %g %g %g %d %g\n",
                    PhydroBody(i)->GlobalID,r,Phydro[i]->Rho,v_r,Phydro[i]->U,
                    SmoothedValue,Weight,
                    Alpha,DAlpha,Phydro[i]->Nlist,SmoothedNumber);

		}
    }
    fclose(fp);
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"cat %s.point.%03d.??? | sort -n > %s.point.%03d",
                name,MPIGetNumProcs(),name,MPIGetNumProcs());
        fprintf(stderr,"%s\n",command);
        system(command);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"rm %s.point.%03d.???",name,MPIGetNumProcs());
        fprintf(stderr,"%s\n",command);
        system(command);
    }



    sprintf(fname,"%s.plot.%03d.%03d",name,MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        double Pos[3] = {PeriodicDistance(PhydroPos(i)[0],0,0),
                         PeriodicDistance(PhydroPos(i)[1],0,1),
                         PeriodicDistance(PhydroPos(i)[2],0,2)};
        r = NORM(Pos);
        fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                r,NORM(PhydroVel(i)),
                    Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->UPred,NORM(Phydro[i]->HydroAcc));
    }
    fclose(fp);

    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"cat %s.plot.%03d.??? | sort -n > %s.plot.%03d",
                name,MPIGetNumProcs(),name,MPIGetNumProcs());
        fprintf(stderr,"%s\n",command);
        system(command);
    }
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        char command[MaxCharactersInLine];
        Snprintf(command,"rm %s.plot.%03d.???",name,MPIGetNumProcs());
        fprintf(stderr,"%s\n",command);
        system(command);
    }

    sprintf(fname,"%s.radial.%03d.%03d",name,MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,fname,"w");
	dr = 0.005;
    r = 0.e0;
	//for(r=dr;r<1.e0;r+=dr){
    int count = 0;
    while(count != Pall.Nhydro){
		rho_r = v_r = u_r = p_r = cs_r = 0.e0;
        double SmoothedValue_r = 0.e0;
        double Weight_r = 0.e0;
		n = 0;
        for(int i=0;i<Pall.Nhydro;i++){
			//r2 = NORM2(PhydroPos(i));
            double Pos[3] = {PeriodicDistance(PhydroPos(i)[0],0,0),
                             PeriodicDistance(PhydroPos(i)[1],0,1),
                             PeriodicDistance(PhydroPos(i)[2],0,2)};
            r2 = NORM2(Pos);
			if( (SQ(r)<=r2) && (r2<SQ(r+dr)) ){
				n++;
				rho_r += Phydro[i]->Rho;
				p_r += Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U;
				cs_r += sqrt(Pall.GGm1*Phydro[i]->U);
				v_r += DOT_PRODUCT(PhydroPos(i),PhydroVel(i))/sqrt(r2);
				u_r += Phydro[i]->U;
#ifdef USE_DISPH //{
                SmoothedValue_r += Phydro[i]->EnergyDensity;
                Weight_r += Phydro[i]->Mass*Phydro[i]->U;
#elif defined(USE_SPSPH) //}//{
                SmoothedValue_r += Phydro[i]->PseudoDensity;
                Weight_r += Phydro[i]->Zw;
#else //}//{
                SmoothedValue_r += Phydro[i]->Rho;
                Weight_r += Phydro[i]->Mass;
#endif // USE_DISPH //}
			}
		}
		if(n != 0){
			fprintf(fp,"%g %g %g %g %g %g %g %g\n",
                r+0.5*dr,rho_r/(double)n,p_r/(double)n,
                    v_r/(double)n,u_r/(double)n,fabs(v_r/cs_r),
                    SmoothedValue_r/(double)n,Weight_r/(double)n);
            count += n;
		}
        r += dr;
	}
	fclose(fp);

	return;
}


void FileIO3DCollapse(void){

    static double IOTiming[] = {0.0,0.8,1.3,1.7,2.6,3.0};
    int Num = sizeof(IOTiming)/sizeof(double);
    static int IOstep = 0;
    char fname[MaxCharactersInLine];

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(IOstep < Num){
            if(IOTiming[IOstep] < Pall.TCurrent){
                sprintf(fname,"io.%02d",IOstep);
                OutPut3DCollapseTest(fname);
                IOstep ++;
            }
        }
    }

    return;
}


void OutPutZeroPlane(char *name){

    //const static double ZeroPlane = 1.0/128.0; 
    const static double ZeroPlane = 1.0/32.0; 
    char fname[MaxCharactersInLine];
    FILE *fp;

    sprintf(fname,"%s.zero",name);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        //if(fabs(PhydroPos(i)[2]) < ZeroPlane){
        if(fabs(PhydroPos(i)[2]) < ZeroPlane){
            fprintf(fp,"%g %g %g %g %g %g %g %g %g %g\n",
                    PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                    PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                    Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->UPred,
                    NORM(Phydro[i]->HydroAcc));
        }
    }
    fclose(fp);

	return;
}

#ifdef TASK_BLAST_WAVE
void FileIOBlastWave(void){

    static double IOTiming[] = {0.0,0.005,0.02,0.04,0.08};
    int Num = sizeof(IOTiming)/sizeof(double);
    static int IOstep = 0;
    char fname[MaxCharactersInLine];

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(IOstep < Num){
            if(IOTiming[IOstep] < Pall.TCurrent){
                sprintf(fname,"BlastWave.%02d",IOstep);
                OutPut3DCollapseTest(fname);
                OutPutZeroPlane(fname);
                OutPutBlastWaveAll(fname);
                IOstep ++;
            }
        }
    }
    return;
}

void OutPutBlastWaveAll(char *name){

    char fname[MaxCharactersInLine];
    FILE *fp;

    // presearch peaks.
    int index = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        //if(Phydro[index]->RhoPred > Phydro[i]->RhoPred)
        if(Phydro[index]->U < Phydro[i]->U)
            index = i;
    }

    sprintf(fname,"%s.all",name);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g %g\n",
                PhydroBody(i)->GlobalID,PhydroPosP(i)[0],PhydroPosP(i)[1],PhydroPosP(i)[2],
                /*
                PhydroBody(i)->GlobalID,
                (PhydroPosP(i)[0]-PhydroPosP(index)[0]),
                (PhydroPosP(i)[1]-PhydroPosP(index)[1]),
                (PhydroPosP(i)[2]-PhydroPosP(index)[2]),
                */
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->RhoPred,Phydro[i]->KernelPred,Phydro[i]->UPred,
                NORM(Phydro[i]->HydroAcc));
    }
    fclose(fp);

    sprintf(fname,"%s.radialpoint",name);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
		double r = DISTANCE(PhydroPosP(i),PhydroPosP(index));
		if( r > 0 ){
			double v_r = DOT_PRODUCT(PhydroPosP(i),PhydroVel(i))/r;
			fprintf(fp,"%ld %e %e %e %e\n",PhydroBody(i)->GlobalID,r,Phydro[i]->RhoPred,v_r,Phydro[i]->UPred);
		} else {
			fprintf(fp,"%ld %e %e %e %e\n",PhydroBody(i)->GlobalID,0.0,Phydro[i]->RhoPred,0.0,Phydro[i]->UPred);
        }
    }
    fclose(fp);

	return;
}
#endif

void FileIOSelfSimilar(void){

    static double IOTiming[] = {0.0,14.0,30.0,60.0,120.0};
    int Num = sizeof(IOTiming)/sizeof(double);
    static int IOstep = 0;
    char fname[MaxCharactersInLine];

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(IOstep < Num){
            if(IOTiming[IOstep] < Pall.TCurrent){
                sprintf(fname,"io.%02d",IOstep);
                OutPut3DCollapseTest(fname);
                IOstep ++;
            }
        }
    }

    return;
}

#ifdef TASK_1D_SHOCKE_TUBE
// Sod Shock Tube Test
void InitShockTube(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    int NBlock = Number >> 2;
    int Ntotal = NBlock << 2;
    fprintf(stderr,"Ntotal,NBlock,Number = %d %d %d\n",Ntotal,NBlock,Number);

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Ntotal;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    dprintlmpi(AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    dprintlmpi(Ntotal);
    int mycount = 0;
    for(int i=0;i<Ntotal;i++){
        if(i%NProcs == MyID){
            PhydroBody(mycount)->Active = ON;
            PhydroBody(mycount)->Use = ON;
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->GlobalID = i;
            if(i<3*NBlock){
                PhydroPos(mycount)[0] = -1.e0+1.0/(3.e0*NBlock)*i;
                Phydro[mycount]->Kernel = 1.e0/((double)NBlock);
                Phydro[mycount]->Rho = 3.e0;
                Phydro[mycount]->U = 3.e0/(Pall.Gm1*Phydro[mycount]->Rho);
            }else{
                PhydroPos(mycount)[0] = 1.e0/(NBlock)*(i-3*NBlock);
                Phydro[mycount]->Kernel = 1.e0/((double)NBlock);
                Phydro[mycount]->Rho = 1.e0;
                Phydro[mycount]->U = 1.e0/(Pall.Gm1*Phydro[mycount]->Rho);
            }
            //PhydroMass(mycount) = 1.e0/((double)Ntotal);
            PhydroBody(mycount)->Mass = PhydroMass(mycount) = 1.e0/((double)NBlock);
            PhydroPos(mycount)[1] = PhydroPos(mycount)[2] = 0.e0;
            PhydroVel(mycount)[0] = PhydroVel(mycount)[1] = PhydroVel(mycount)[2] = 0.e0;
            Phydro[mycount]->HydroAcc[0] = Phydro[mycount]->HydroAcc[1] = Phydro[mycount]->HydroAcc[2] = 0.e0;
            Phydro[mycount]->Use = ON; 
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->Eps = 1.0;
            mycount ++;
        }
    }
    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = Ntotal;

    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 5;
    Pall.Npm = 1;

    Pall.TEnd = 0.1;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/ST1D.ASCII");
    strcpy(Pall.BaseFileName,"./data/ST1D");
    strcpy(Pall.RestartFileName,"./data/ST1D.dump");

    // InitLogFiles();
    // char fname[MaxCharactersInLine];
    // sprintf(fname,"InitShockTube.%02d.data",MyID);
    // OutPutShockTube(fname);

    return;
}

// Sod shock tube 
// Ref. Sod (1978)
void InitSodShockTube(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    int NBlock = Number/5;
    int NLeft = 4*NBlock;
    int NRight = NBlock;
    int Ntotal = NLeft+NRight;

    fprintf(stderr,"Ntotal, NBlock = %d, %d\n",Ntotal,NBlock);

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Ntotal;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    dprintlmpi(AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else // USE_VARIABLE_ALPHA //}//{
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    double dx_l = 0.5/((double)NLeft);
    double dx_r = 0.5/((double)NRight);

    dprintlmpi(Ntotal);
    int mycount = 0;
    for(int i=0;i<Ntotal;i++){
        if(i%NProcs == MyID){
            PhydroBody(mycount)->Active = ON;
            PhydroBody(mycount)->Use = ON;
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->GlobalID = i;
            if(i<NLeft){
                PhydroPos(mycount)[0] = -5.e-1+i*dx_l;
                Phydro[mycount]->Kernel = 2.0*dx_l;
                Phydro[mycount]->Rho = 1.e0;
                Phydro[mycount]->U = 1.e0/(Pall.Gm1*Phydro[mycount]->Rho);
            }else{
                PhydroPos(mycount)[0] = (i-NLeft)*dx_r;
                Phydro[mycount]->Kernel = 2.0*dx_r;
                Phydro[mycount]->Rho = 0.25;
                Phydro[mycount]->U = 0.1795/(Pall.Gm1*Phydro[mycount]->Rho);
            }
#ifdef EVALUATE_KERNEL_BY_ITERATION
            Phydro[mycount]->Kernel = 1.e0/((double)NBlock);
#endif
            PhydroBody(mycount)->Mass = PhydroMass(mycount) = 0.5*1.25/((double)Ntotal);
            PhydroPos(mycount)[1] = PhydroPos(mycount)[2] = 0.e0;
            PhydroVel(mycount)[0] = PhydroVel(mycount)[1] = PhydroVel(mycount)[2] = 0.e0;
            Phydro[mycount]->HydroAcc[0] = Phydro[mycount]->HydroAcc[1] = Phydro[mycount]->HydroAcc[2] = 0.e0;
            Phydro[mycount]->Use = ON; 
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->Eps = 1.0;
            mycount ++;
        }
    }
    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = Ntotal;


    ActivateAllparticles();


    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 6;
    Pall.Npm = 2;

    Pall.TEnd = 0.1;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/ST1D.ASCII");
    strcpy(Pall.BaseFileName,"./data/ST1D");
    strcpy(Pall.RestartFileName,"./data/ST1D.dump");

    // InitLogFiles();
    // char fname[MaxCharactersInLine];
    // sprintf(fname,"InitShockTube.%02d.data",MyID);
    // OutPutShockTube(fname);

    // for(int i=0;i<Pall.Nhydro;i++){
        // Phydro[i]->Kernel = Phydro[i]->KernelPred = gsl_rng_uniform(RandomGenerator)*Phydro[i]->Kernel;
    // }

    return;
}

// Sod Shock Tube Test (Hernquist Katz '89)
void InitShockTubeHK(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int NBlock = Number / 5;
    int Ntotal = NBlock * 5;
    fprintf(stderr,"Ntotal,NBlock,Number = %d %d %d\n",Ntotal,NBlock,Number);

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitLength = Pall.UnitTime = Pall.UnitMass = Pall.GravConst = 1.0;

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Ntotal;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    // Allocate Particle Data Structures.

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    int mycount = 0;
    for(int i=0;i<Ntotal;i++){
        if(i%NProcs == MyID){
            PhydroBody(mycount)->Active = ON;
            PhydroBody(mycount)->Use = ON;
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->GlobalID = i;
            if(i<4*NBlock){
                PhydroPos(mycount)[0] = -1.e0+(1.0/(4.e0*NBlock))*i;
                Phydro[mycount]->Kernel = 1.e0/((double)NBlock);
                Phydro[mycount]->Rho = 1.e0;
                Phydro[mycount]->U = 1.0/(Pall.Gm1*Phydro[mycount]->Rho);
            }else{
                PhydroPos(mycount)[0] = 1.e0/(NBlock)*(i-4*NBlock);
                Phydro[mycount]->Kernel = 1.e0/((double)NBlock);
                Phydro[mycount]->Rho = 0.25;
                Phydro[mycount]->U = 0.1795/(Pall.Gm1*Phydro[mycount]->Rho);
            }
            //PhydroMass(mycount) = 0.75/((double)Ntotal);
            PhydroMass(mycount) = 1.25/((double)Ntotal);
            PhydroPos(mycount)[1] = PhydroPos(mycount)[2] = 0.e0;
            PhydroVel(mycount)[0] = PhydroVel(mycount)[1] = PhydroVel(mycount)[2] = 0.e0;
            Phydro[mycount]->HydroAcc[0] = Phydro[mycount]->HydroAcc[1] = Phydro[mycount]->HydroAcc[2] = 0.e0;
            Phydro[mycount]->Use = ON; 
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->Eps = 0.1;
#ifdef USE_VARIABLE_ALPHA
            //Phydro[mycount]->Alpha = 0.01;
            double h2 = 0.005;
            Phydro[i]->Alpha = Pall.ViscousAlphaMax*exp(-SQ(PhydroPos(i)[0])/h2);
            if(Phydro[i]->Alpha < Pall.ViscousAlphaMin)
                Phydro[i]->Alpha = Pall.ViscousAlphaMin;
#endif // USE_VARIABLE_ALPHA
            mycount ++;
        }
    }
    Pall.Ntotal = Pall.Nhydro = Pall.NActives = Pall.NActivesHydro =  mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = Pall.NActives_t = Pall.NActivesHydro_t = Ntotal;

    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 5;
    //Pall.Npm = 1;
    Pall.Npm = 1;

    Pall.TEnd = 0.15;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/ST.ASCII");
    strcpy(Pall.BaseFileName,"./data/ST");
    strcpy(Pall.RestartFileName,"./data/ST.dump");

    return;
}

void InitShockTube123Problem(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int NBlock = Number / 2;
    int Ntotal = NBlock * 2;
    fprintf(stderr,"Ntotal,NBlock,Number = %d %d %d\n",Ntotal,NBlock,Number);

    memset(&Pall,0,sizeof(struct StructPall));

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Ntotal;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    dprintlmpi(AllocationSize);
    
    PbodySize = AllocationSize;
    PhydroSize = AllocationSize;
    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(AllocationSize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(AllocationSize*sizeof(StructPhydroptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    // Allocate Particle Data Structures.

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    int mycount = 0;
    for(int i=0;i<Ntotal;i++){
        if(i%NProcs == MyID){
            PhydroBody(mycount)->Active = ON;
            PhydroBody(mycount)->Use = ON;
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->GlobalID = i;
            if(i<NBlock){
                PhydroPos(mycount)[0] = -1.e0+(1.0/NBlock)*i;
                Phydro[mycount]->Kernel = 1.e0/((double)NBlock);
                Phydro[mycount]->Rho = 1.e0;
                Phydro[mycount]->U = 0.4/(Pall.Gm1*Phydro[mycount]->Rho);
                PhydroVel(mycount)[0] = -2.0;
            }else{
                PhydroPos(mycount)[0] = 1.e0/(NBlock)*(i-NBlock);
                Phydro[mycount]->Kernel = 1.e0/((double)NBlock);
                Phydro[mycount]->Rho = 1.e0;
                Phydro[mycount]->U = 0.4/(Pall.Gm1*Phydro[mycount]->Rho);
                PhydroVel(mycount)[0] = 2.0;
            }
            //PhydroMass(mycount) = 0.75/((double)Ntotal);
            PhydroMass(mycount) = 2.0/((double)Ntotal);
            PhydroPos(mycount)[1] = PhydroPos(mycount)[2] = 0.e0;
            PhydroVel(mycount)[1] = PhydroVel(mycount)[2] = 0.e0;
            Phydro[mycount]->HydroAcc[0] = Phydro[mycount]->HydroAcc[1] = Phydro[mycount]->HydroAcc[2] = 0.e0;
            Phydro[mycount]->Use = ON; 
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->Eps = 0.1;
            mycount ++;
        }
    }
    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = Ntotal;

    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 5;
    Pall.Npm = 1;

    Pall.TEnd = 0.15;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/ST1D.ASCII");
    strcpy(Pall.BaseFileName,"./data/ST1D");
    strcpy(Pall.RestartFileName,"./data/ST1D.dump");

    return;
}

void OutPutShockTube(char *fname){

    FILE *fp;
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
#ifdef USE_VARIABLE_ALPHA //{
        double Alpha = Phydro[i]->Alpha;
#else // USE_VARIABLE_ALPHA //}
        double Alpha = Pall.HydroAlpha;
#endif // USE_VARIABLE_ALPHA}

#ifdef USE_DISPH //{
        double Fund = Phydro[i]->EnergyDensity;
        double P = Pall.Gm1*Phydro[i]->EnergyDensity;
        double U = Phydro[i]->U;
        double Entropy = Pall.Gm1*Phydro[i]->U*pow(Phydro[i]->Rho,1-Pall.Gamma);
#elif defined(USE_SPSPH) //{
        double Fund = Phydro[i]->PseudoDensity;
        double P = Pall.Gm1*Phydro[i]->Mass*Phydro[i]->U*Phydro[i]->PseudoDensity/(Phydro[i]->Zw);
        double U = Phydro[i]->U;
        double Entropy = Pall.Gm1*Phydro[i]->U*pow(Phydro[i]->Rho,1-Pall.Gamma);
#else  //USE_DISPH //}//{
        double Fund = Phydro[i]->Rho;
        double P = Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U;
        double U = Phydro[i]->U;
        double Entropy = Pall.Gm1*Phydro[i]->U*pow(Phydro[i]->Rho,1-Pall.Gamma);
#endif  //USE_DISPH //}

        fprintf(fp,"%ld %d %e %e %e %e %e %e %e %e %e %e\n",PhydroBody(i)->GlobalID,Phydro[i]->Nlist,
                PhydroPos(i)[0],PhydroVel(i)[0],Phydro[i]->HydroAcc[0],Phydro[i]->Rho,Phydro[i]->Kernel,Fund,U,Entropy,P,Alpha);
#if 0
#ifdef USE_DISPH
#ifdef USE_VARIABLE_ALPHA
        fprintf(fp,"%ld %d %e %e %e %e %e %e %e %e %e\n",PhydroBody(i)->GlobalID,Phydro[i]->Nlist,
                PhydroPos(i)[0],Phydro[i]->EnergyDensity/Phydro[i]->U,Phydro[i]->Kernel,
            Pall.Gm1*Phydro[i]->EnergyDensity,Phydro[i]->U,PhydroVel(i)[0],Phydro[i]->HydroAcc[0],
                Pall.Gm1*Phydro[i]->U/pow(Phydro[i]->Rho,Pall.Gm1),Phydro[i]->Alpha);
#else // USE_VARIABLE_ALPHA
        fprintf(fp,"%ld %d %e %e %e %e %e %e %e %e\n",PhydroBody(i)->GlobalID,Phydro[i]->Nlist,
                PhydroPos(i)[0],Phydro[i]->EnergyDensity/Phydro[i]->U,Phydro[i]->Kernel,
            Pall.Gm1*Phydro[i]->EnergyDensity,Phydro[i]->U,PhydroVel(i)[0],Phydro[i]->HydroAcc[0],
                Pall.Gm1*Phydro[i]->U/pow(Phydro[i]->Rho,Pall.Gm1));
#endif // USE_VARIABLE_ALPHA
#else // USE_DISPH
#ifdef USE_VARIABLE_ALPHA
        fprintf(fp,"%ld %d %e %e %e %e %e %e %e %e %e\n",PhydroBody(i)->GlobalID,Phydro[i]->Nlist,
                PhydroPos(i)[0],Phydro[i]->Rho,Phydro[i]->Kernel,
            Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U,Phydro[i]->U,PhydroVel(i)[0],Phydro[i]->HydroAcc[0],
                Pall.Gm1*Phydro[i]->U/pow(Phydro[i]->Rho,Pall.Gm1),Phydro[i]->Alpha);
#else // USE_VARIABLE_ALPHA
        fprintf(fp,"%ld %d %e %e %e %e %e %e %e %e\n",PhydroBody(i)->GlobalID,Phydro[i]->Nlist,
                PhydroPos(i)[0],Phydro[i]->Rho,Phydro[i]->Kernel,
            Pall.Gm1*Phydro[i]->Rho*Phydro[i]->U,Phydro[i]->U,PhydroVel(i)[0],Phydro[i]->HydroAcc[0],
                Pall.Gm1*Phydro[i]->U/pow(Phydro[i]->Rho,Pall.Gm1));
#endif // USE_VARIABLE_ALPHA
#endif // USE_DISPH
#endif
    }
    fclose(fp);

	return;
}
#endif //TASK_1D_SHOCKE_TUBE

#ifdef TASK_1D_TWOFLUIDS //{
void Init1DTwoFluids(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else // USE_VARIABLE_ALPHA //}//{
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    double mu1 = 1.0;
    double mu2 = 2.0;
    double T = 1.0; double P = 1.0;

    double dx = 1.0/((double)Number);
    double mass = 1.0/((double)Number);

    for(int i=0;i<Number;i++){
        PhydroBody(i)->Active = ON;
        PhydroBody(i)->Use = ON;
        PhydroBody(i)->Type = TypeHydro;
        PhydroBody(i)->GlobalID = i;

        PhydroPos(i)[0] = dx*i;
        if( (PhydroPos(i)[0]>0.25) && (PhydroPos(i)[0]<0.75) ){
            Phydro[i]->Kernel = 2.0*dx;
            Phydro[i]->Rho = 1.e0*mu1;
            PhydroBody(i)->Mass = mass*mu1;
            Phydro[i]->Tag = 0;
        } else {
            Phydro[i]->Kernel = 2.0*dx;
            Phydro[i]->Rho = 1.e0*mu2;
            PhydroBody(i)->Mass = mass*mu2;
            Phydro[i]->Tag = 1;
        }
        Phydro[i]->U = P/(Pall.Gm1*Phydro[i]->Rho);

        PhydroPos(i)[1] = PhydroPos(i)[2] = 0.e0;
        PhydroVel(i)[0] = PhydroVel(i)[1] = PhydroVel(i)[2] = 0.e0;
        PhydroAcc(i)[0] = PhydroAcc(i)[1] = PhydroAcc(i)[2] = 0.e0;
        Phydro[i]->HydroAcc[0] = Phydro[i]->HydroAcc[1] = Phydro[i]->HydroAcc[2] = 0.e0;
        Phydro[i]->Use = ON; 
        PhydroBody(i)->Type = TypeHydro;
        PhydroBody(i)->Eps = 1.0;
        //fprintf(stderr,"%g %g %g\n",PhydroPos(i)[0],Phydro[i]->Rho, );
    }


    Pall.Ntotal = Pall.Nhydro = Number;
    Pall.Ntotal_t = Pall.Nhydro_t = Number;


    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 10;
    Pall.Npm = 4;

    Pall.TEnd = 100.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.01*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/TF1D.ASCII");
    strcpy(Pall.BaseFileName,"./data/TF1D");
    strcpy(Pall.RestartFileName,"./data/TF1D.dump");

    return;
}
#endif // TASK_1D_TWOFLUIDS //}

#ifdef TASK_KEPLER //{
#define GM (1.0)
void InitKepler(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);


    int p = 6;
    int Nring = Number;
    double r_in = 0.5;
    double r_out = 2;
    double dr = r_out/Nring;
    int i_min = (r_in/dr)+1;

    int count = 0; int mycount = 0;
    for(int i=i_min;i<Nring;i++){
        double r = i*dr;
        double phi = 2*M_PI/(i*p);
        for(int k=0;k<i*p;k++){
            if(count%NProcs == MyID)
                mycount ++;
            count ++;
        }
    }
    dprintlmpi(count);
    dprintlmpi(mycount);


    // Allocate Particle Data Structures.
    int AllocationSize = mycount;
    // for(int i=0;i<Number;i++)
        // if(i%NProcs == MyID)
            // AllocationSize ++;

    //Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    //Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else // USE_VARIABLE_ALPHA //}//{
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    //double Surface = 4.0*4.0;
    double Surface = M_PI*(2*2-0.5*0.5);
    double Rho = 1.0;
    double TotalMass = Surface*Rho;
    double Mass = TotalMass/((double)count);
    double Kernel = 10*sqrt(Surface)/count;

    double P = 1.e-6;
    double u = P/(Pall.Gm1*Rho);

    count = 0;
    mycount = 0;
    for(int i=i_min;i<Nring;i++){

        double r = i*dr;
        double phi = 2*M_PI/(i*p);
        for(int k=0;k<i*p;k++){

            double x = r*cos(phi*k);
            double y = r*sin(phi*k);

            double r = sqrt(x*x+y*y);
            double Vc = sqrt(GM/r);

            if(count%NProcs == MyID){

                PhydroBody(mycount)->Active = ON;
                PhydroBody(mycount)->Use = ON;
                PhydroBody(mycount)->Type = TypeHydro;
                PhydroBody(mycount)->GlobalID = count;


                Phydro[mycount]->Kernel = Kernel;
                Phydro[mycount]->Rho = 1.e0;
                PhydroBody(mycount)->Mass = Mass;
                Phydro[mycount]->U = u;
                Phydro[mycount]->Tag = 0;

                PhydroPos(mycount)[0] = x;
                PhydroPos(mycount)[1] = y;
                PhydroVel(mycount)[0] = -Vc*(y/r);
                PhydroVel(mycount)[1] = +Vc*(x/r);

                PhydroPos(mycount)[2] = 0.e0;
                PhydroVel(mycount)[2] = 0.e0;
                PhydroAcc(mycount)[0] = PhydroAcc(mycount)[1] = PhydroAcc(mycount)[2] = 0.e0;
                Phydro[mycount]->HydroAcc[0] = Phydro[mycount]->HydroAcc[1] = Phydro[mycount]->HydroAcc[2] = 0.e0;
                Phydro[mycount]->Use = ON; 
                PhydroBody(mycount)->Type = TypeHydro;
                PhydroBody(mycount)->Eps = 1.0;


                mycount ++;
            }
            count ++;
        }
    }

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;


    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 64;
    Pall.Npm = 4;

    Pall.TEnd = 600.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 1;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/kp.ASCII");
    strcpy(Pall.BaseFileName,"./data/kp");
    strcpy(Pall.RestartFileName,"./data/kp.dump");

    return;
}


void KeplerPotential(void){

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            double r = NORM(Pbody[i]->PosP);
            if(r < 0.25){
                r = sqrt(SQ(r) + SQ(0.25));
                Pbody[i]->Acc[0] += -(GM/CUBE(r))*Pbody[i]->PosP[0];
                Pbody[i]->Acc[1] += -(GM/CUBE(r))*Pbody[i]->PosP[1];
            } else {
                Pbody[i]->Acc[0] += -(GM/CUBE(r))*Pbody[i]->PosP[0];
                Pbody[i]->Acc[1] += -(GM/CUBE(r))*Pbody[i]->PosP[1];
            }
        }
    }

    return ;
}
#endif // TASK_KEPLER //}

// Cold Collapse Test.

/*
 * This function returns a double precision value following the Gaussian with
 * the mean value of 0.0 and dispersion = 1.0 by Box-Muller method. Before the
 * use (calling) of this function, you must initialize the GSL random
 * generator.
 */
static double Gaussian(void){

    double x,y,r2;

    do{
        x = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        y = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        r2 = SQ(x)+SQ(y);
    } while ((r2 >= 1.0) || (r2 == 0.0));

    return (sqrt(-2.0*log(r2)/r2)*x);
}

#ifdef TASK_COLD_COLLAPSE
static void CalcPotentialDirectForColdCollaspe(void);

/*
 * This function needs a total number of particles and a virial ratio in order
 * to decide velocity dispersion.
 */ 
void InitColdCollapseTest(const int Number, const double rv, const int MakePotential){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);
    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1;
    //Pall.GravConst = GetUnitGravitationalConstant();
    Pall.GravConst = 1;

    // Allocate Particle Data Structures.
    int mycount = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            mycount ++;

    Pall.Ntotal_t = Number;
    Pall.Ntotal = mycount;

    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    // Allocate Particle Data Structures.

    // Generate a particles distribution.
    mycount = 0; double Pos[3];
    double eps = 1.0/(sqrt((double)Number));
    for(int i=0;i<Number;i++){
        if(i%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeDM;
            Pbody[mycount]->GlobalID = i;

            Pbody[mycount]->Mass = 1.e0/((double)Number);
            Pbody[mycount]->Eps = eps;

            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
            Pbody[mycount]->Pos[0] = Pos[0];
            Pbody[mycount]->Pos[1] = Pos[1];
            Pbody[mycount]->Pos[2] = Pos[2];

            Pbody[mycount]->InteractionList = 1;

            mycount ++;
        } else {
            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
        }
    }
    dprintlmpi(mycount);


    double sigma;
    if(MakePotential == ON){
        //PotentialDirect();
        CalcPotentialDirectForColdCollaspe();

        // Calc Potential.
        double W = 0.e0;
        for(int i=0;i<mycount;i++)
            W += Pbody[i]->Pot;
        eprintlmpi(W);
        double GW;
        MPI_Allreduce(&W,&GW,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        W = GW;

        double Mtotal = 1.e0;
        sigma = sqrt(2.0*fabs(W)*rv/(3.0*Mtotal));

        eprintlmpi(sigma);
    } else {
        sigma = 0.e0;
    }

    // Calc sigma and add velocity for particles.
    mycount = 0;
    double Vel[3];
    for(int i=0;i<Number;i++){
        if(i%NProcs == MyID){
            Pbody[mycount]->Vel[0] = sigma*Gaussian();
            Pbody[mycount]->Vel[1] = sigma*Gaussian();
            Pbody[mycount]->Vel[2] = sigma*Gaussian();
            mycount ++;
        } else {
            Vel[0] = gsl_ran_gaussian(RandomGenerator,sigma);
            Vel[1] = gsl_ran_gaussian(RandomGenerator,sigma);
            Vel[2] = gsl_ran_gaussian(RandomGenerator,sigma);
        }
    }

    /*
    for(int i=0;i<Number;i++){
        fprintf(stderr,"%d %g %g %g %g %g %g\n",i,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
    }
    */

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 1.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 10*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/CC.ASCII");
    strcpy(Pall.BaseFileName,"./data/CC");
    strcpy(Pall.RestartFileName,"./data/CC.dump");

    return;
}

/*
 * This function makes an initial condition for a cold collapse test with two
 * different mass particles. A half of particles have the mass of m1, and others
 * have the mass of m2. Since the total mass of the system is 1, the particle
 * masses should satisfy the relation: 1 = N(m1+m2)/2. When we adopt m1:m2 =
 * 1:2, the relation is 1 = 3 N m1/2. Thus m1 = 2/(3N). The third parameter of
 * this function gives m1/m2.
 */
void InitColdCollapseTestMixed(const int Number, const double rv, const double fraction, const int MakePotential){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);
    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1;
    //Pall.GravConst = GetUnitGravitationalConstant();
    Pall.GravConst = 1;

    // Allocate Particle Data Structures.
    int mycount = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            mycount ++;

    Pall.Ntotal_t = Number;
    Pall.Ntotal = mycount;

    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    // Allocate Particle Data Structures.
    
    
    double m1 = 2.0*fraction/((double)Number*(fraction+1.0));
    double m2 = 2.0/((double)Number*(fraction+1.0));
    fprintf(stderr,"Mass 1 and 2 are %g and %g\n",m1,m2); 

    double eps = 1.0/(sqrt((double)Number));
    double m0 = 1.0/(double)Number;
    double eps1 = eps*cbrt(m1/m0);
    double eps2 = eps*cbrt(m2/m0);

    // Generate a particles distribution.
    mycount = 0;


    for(int i=0;i<Number;i++){
        double Pos[3];
        if(i%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeDM;
            Pbody[mycount]->GlobalID = i;

            if(i%2 == 0){
                Pbody[mycount]->Mass = m1;
                Pbody[mycount]->Eps = eps1;
            }else{
                Pbody[mycount]->Mass = m2;
                Pbody[mycount]->Eps = eps2;
            }

            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
            Pbody[mycount]->Pos[0] = Pos[0];
            Pbody[mycount]->Pos[1] = Pos[1];
            Pbody[mycount]->Pos[2] = Pos[2];

            Pbody[mycount]->InteractionList = 1;

            mycount ++;
        } else {
            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
        }
    }
    dprintlmpi(mycount);

    double sigma;
    if(MakePotential == ON){
        //PotentialDirect();
        CalcPotentialDirectForColdCollaspe();

        // Calc Potential.
        double W = 0.e0;
        for(int i=0;i<mycount;i++)
            W += Pbody[i]->Pot;
        eprintlmpi(W);
        double GW;
        MPI_Allreduce(&W,&GW,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        W = GW;

        double Mtotal = 1.e0;
        sigma = sqrt(2.0*fabs(W)*rv/(3.0*Mtotal));

        eprintlmpi(sigma);
    } else {
        sigma = 0.e0;
    }

    // Calc sigma and add velocity for particles.
    mycount = 0;
    double Vel[3];
    for(int i=0;i<Pall.Ntotal;i++){
        if(i%NProcs == MyID){
            Pbody[mycount]->Vel[0] = sigma*Gaussian();
            Pbody[mycount]->Vel[1] = sigma*Gaussian();
            Pbody[mycount]->Vel[2] = sigma*Gaussian();
            mycount ++;
        } else {
            Vel[0] = gsl_ran_gaussian(RandomGenerator,sigma);
            Vel[1] = gsl_ran_gaussian(RandomGenerator,sigma);
            Vel[2] = gsl_ran_gaussian(RandomGenerator,sigma);
        }
    }

    double total_mass = 0.e0;
    double MeanVel[3] = {0.e0,0.e0,0.e0};
    for(int i=0;i<Pall.Ntotal;i++){
        MeanVel[0] += Pbody[i]->Mass*Pbody[i]->Vel[0];
        MeanVel[1] += Pbody[i]->Mass*Pbody[i]->Vel[1];
        MeanVel[2] += Pbody[i]->Mass*Pbody[i]->Vel[2];
        total_mass += Pbody[i]->Mass;
    }
    //gprintlmpi(total_mass);
    //fflush(NULL);
    MeanVel[0] /= total_mass;
    MeanVel[1] /= total_mass;
    MeanVel[2] /= total_mass;
    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Vel[0] -= MeanVel[0];
        Pbody[i]->Vel[1] -= MeanVel[1];
        Pbody[i]->Vel[2] -= MeanVel[2];
    }

    /*
    for(int i=0;i<Number;i++){
        fprintf(stderr,"%d %g %g %g %g %g %g %g %g\n",i,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Mass,Pbody[i]->Eps);
    }
    exit(1);
    */

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 1.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 10*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/CC.ASCII");
    strcpy(Pall.BaseFileName,"./data/CC");
    strcpy(Pall.RestartFileName,"./data/CC.dump");

    return;
}

static void CalcPotentialDirectForColdCollaspe(void){

    int NumProcs = MPIGetNumProcs();
    double r,r2,xij[3],EPS;
    MPI_Status  mpi_status;

    /* Self-Domain Potential */
    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Pot = 0.e0;
        for(int j=0;j<Pall.Ntotal;j++){
            xij[0] = Pbody[i]->Pos[0]-Pbody[j]->Pos[0];
            xij[1] = Pbody[i]->Pos[1]-Pbody[j]->Pos[1];
            xij[2] = Pbody[i]->Pos[2]-Pbody[j]->Pos[2];
            r2 = NORM2(xij);
            EPS = 0.5*(Pbody[i]->Eps+Pbody[j]->Eps);
            //if(r2 < SQ(TINY*EPS)) continue;

            r = sqrt(r2+SQ(EPS));
            Pbody[i]->Pot -= Pbody[j]->Mass/r;
        }
    }

    struct StructPotential{
        double Pos[3];
        double Eps;
        double Mass;
    } *PotentialExportSend,*PotentialExportRecv;

    /* Get Other Domains and Calculate Gravity */
    for(int i=0;i<NumProcs-1;i++){
        int SendThisTime,RecvThisTime;
        int SendAlready = 0;
        int RecvAlready = 0;
        int SendData = CommunicationTable[i].SendSize;
        int RecvData = CommunicationTable[i].RecvSize;
        while((SendData != 0)&&(RecvData != 0)){
            SendThisTime = CommunicationTable[i].SendSize-SendAlready;
            CheckSizeofBufferExportSendIndex(SendThisTime,sizeof(struct StructPotential),0);
            PotentialExportSend = BufferExportSend[0];

            int NActives = 0;
            while(NActives != SendThisTime){
                PotentialExportSend[NActives].Pos[0] = Pbody[NActives+SendAlready]->Pos[0];
                PotentialExportSend[NActives].Pos[1] = Pbody[NActives+SendAlready]->Pos[1];
                PotentialExportSend[NActives].Pos[2] = Pbody[NActives+SendAlready]->Pos[2];
                PotentialExportSend[NActives].Mass = Pbody[NActives+SendAlready]->Mass;
                PotentialExportSend[NActives].Eps = Pbody[NActives+SendAlready]->Eps;
                NActives ++;
            }
            RecvThisTime = CommunicationTable[i].RecvSize-RecvAlready;
            CheckSizeofBufferExportRecv(RecvThisTime,sizeof(struct StructPotential));
            PotentialExportRecv = BufferExportRecv;

            // Communication
            MPI_Sendrecv(PotentialExportSend,SendThisTime*sizeof(struct StructPotential),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_POTENTIAL_DIRECT_SENDRECV,
                PotentialExportRecv,RecvThisTime*sizeof(struct StructPotential),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_POTENTIAL_DIRECT_SENDRECV,
                        MPI_COMM_WORLD,&mpi_status);

            // Calculation
            for(int j=0;j<Pall.Ntotal;j++){
                for(int k=0;k<RecvThisTime;k++){
                    xij[0] = Pbody[j]->Pos[0]-PotentialExportRecv[k].Pos[0];
                    xij[1] = Pbody[j]->Pos[1]-PotentialExportRecv[k].Pos[1];
                    xij[2] = Pbody[j]->Pos[2]-PotentialExportRecv[k].Pos[2];
                    r2 = NORM2(xij);
                    EPS = 0.5*(Pbody[j]->Eps+PotentialExportRecv[k].Eps);
                    //if(r2 < SQ(TINY*EPS)) continue;

                    r = sqrt(r2+SQ(EPS));
                    Pbody[j]->Pot -= PotentialExportRecv[k].Mass/r;
                }
            }

            SendData -= SendThisTime;
            RecvData -= RecvThisTime;
            SendAlready += SendThisTime;
            RecvAlready += RecvThisTime;
        }
    }


    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Pot += Pbody[i]->Mass/Pbody[i]->Eps;
        Pbody[i]->Pot *= 0.5*Pall.GravConst*Pbody[i]->Mass;
    }

    if(NumProcs > 1){
        Delete(PotentialExportSend);
        Delete(PotentialExportRecv);
    }

    return;
}
#endif

#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
/*
 * This function makes a spherical distribution of particles where particles
 * with two different mass consist of the sphere. A half of particles have the
 * mass of m1, and others have the mass of m2. Since the total mass of the
 * system is 1, the particle masses should satisfy the relation: 1 = N(m1+m2)/2.
 * When we adopt m1:m2 = 1:2, the relation is 1 = 3 N m1/2. Thus m1 = 2/(3N).
 * The 3rd parameter provides the radius of the sphere, while the 4th parameter
 * gives tha position X where the gravitaional potential is evaluated.
 */
void InitSphericalParticleDistributionWithTwoTypes(const int Number, const double fraction, const double Radius, const double Distance){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);
    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1;
    Pall.GravConst = 1;

    // Allocate Particle Data Structures.
    int mycount = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            mycount ++;
    if(MPIGetMyID()==MPI_ROOT_RANK)
        mycount += 1;

    Pall.Ntotal_t = Number+1;
    Pall.Ntotal = mycount;

    int AllocationSize = 1.2*mycount; 
    GenerateStructPbody(AllocationSize);
    // Allocate Particle Data Structures.
    
    double m1 = 2.0*fraction/((double)Number*(fraction+1.0));
    double m2 = 2.0/((double)Number*(fraction+1.0));
    fprintf(stderr,"Mass 1 and 2 are %g and %g\n",m1,m2); 

    double eps = 1.0/(sqrt((double)1024));
    double m0 = 1.0/(double)1024;
    double eps1 = eps*cbrt(m1/m0);
    double eps2 = eps*cbrt(m2/m0);


    // Generate a particles distribution.
    mycount = 0;


    for(int i=0;i<Number;i++){
        double Pos[3];
        if(i%NProcs == MyID){
            Pbody[mycount]->Active = OFF;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeDM;
            Pbody[mycount]->GlobalID = i;

            if(i%2 == 0){
                Pbody[mycount]->Mass = m1;
                Pbody[mycount]->Eps = eps1;
            }else{
                Pbody[mycount]->Mass = m2;
                Pbody[mycount]->Eps = eps2;
            }

            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
            Pbody[mycount]->Pos[0] = Radius*Pos[0];
            Pbody[mycount]->Pos[1] = Radius*Pos[1];
            Pbody[mycount]->Pos[2] = Radius*Pos[2];

            Pbody[mycount]->Vel[0] = 
            Pbody[mycount]->Vel[1] = 
            Pbody[mycount]->Vel[2] = 0.0;

            Pbody[mycount]->InteractionList = 1;

            mycount ++;
        } else {
            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
        }
    }
    dprintlmpi(mycount);

    ///////////////////////////////////////////////////////////
    if(MPIGetMyID()==MPI_ROOT_RANK){
        Pbody[mycount]->Pos[0] = Distance;
        Pbody[mycount]->Pos[1] = 0.e0;
        Pbody[mycount]->Pos[2] = 0.e0;

        Pbody[mycount]->Active = ON;
        Pbody[mycount]->Use = ON;
        Pbody[mycount]->Type = TypeDM;
        Pbody[mycount]->GlobalID = mycount;

        Pbody[mycount]->Mass = m1;
        Pbody[mycount]->Eps = eps1;
    }


    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 1.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 10*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/CC.ASCII");
    strcpy(Pall.BaseFileName,"./data/CC");
    strcpy(Pall.RestartFileName,"./data/CC.dump");

    return;
}

void RestoreSphericalParticleDistribution(const int Number, double Pos[][3], double Mass[], double Eps[], const double M_p, const double Eps_p, const double Radius, const double Distance){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);
    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1;
    Pall.GravConst = 1;

    // Allocate Particle Data Structures.
    int mycount = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            mycount ++;
    if(MPIGetMyID()==MPI_ROOT_RANK)
        mycount += 1;

    Pall.Ntotal_t = Number+1;
    Pall.Ntotal = mycount;

    int AllocationSize = 1.2*mycount; 
    GenerateStructPbody(AllocationSize);
    // Allocate Particle Data Structures.
    
    // Generate a particles distribution.
    mycount = 0;
    for(int i=0;i<Number;i++){
        Pbody[mycount]->Active = OFF;
        Pbody[mycount]->Use = ON;
        Pbody[mycount]->Type = TypeDM;
        Pbody[mycount]->GlobalID = i;

        Pbody[mycount]->Pos[0] = Pos[i][0];
        Pbody[mycount]->Pos[1] = Pos[i][1];
        Pbody[mycount]->Pos[2] = Pos[i][2];
        Pbody[mycount]->Mass = Mass[i];
        Pbody[mycount]->Eps = Eps[i];

        Pbody[mycount]->Vel[0] = 
        Pbody[mycount]->Vel[1] = 
        Pbody[mycount]->Vel[2] = 0.0;

        Pbody[mycount]->InteractionList = 1;

        mycount ++;
    }
    dprintlmpi(mycount);

    ///////////////////////////////////////////////////////////
    if(MPIGetMyID()==MPI_ROOT_RANK){
        Pbody[mycount]->Pos[0] = Distance;
        Pbody[mycount]->Pos[1] = 0.e0;
        Pbody[mycount]->Pos[2] = 0.e0;

        Pbody[mycount]->Active = ON;
        Pbody[mycount]->Use = ON;
        Pbody[mycount]->Type = TypeDM;
        Pbody[mycount]->GlobalID = mycount;

        Pbody[mycount]->Mass = M_p;
        Pbody[mycount]->Eps = Eps_p;
    }


    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 1.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 10*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/CC.ASCII");
    strcpy(Pall.BaseFileName,"./data/CC");
    strcpy(Pall.RestartFileName,"./data/CC.dump");

    return;
}


/*
 * This function makes a particle distribution with an uniform sphere of the
 * radius ``Radius''. There are two groups with different masses. The mass
 * ratio for particle 1 and particle 2 is given by the paramter ``fraction''.
 * The number of particles is ``Number''.
 */
void InitSphericalParticleDistributionForErrorEstimation(const int Number, const double fraction, const double Radius){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);
    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1;
    Pall.GravConst = 1;

    // Allocate Particle Data Structures.
    int mycount = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            mycount ++;

    Pall.Ntotal_t = Number;
    Pall.Ntotal = mycount;

    int AllocationSize = 1.2*mycount; 
    GenerateStructPbody(AllocationSize);
    // Allocate Particle Data Structures.
    
    double m1 = 2.0*fraction/((double)Number*(fraction+1.0));
    double m2 = 2.0/((double)Number*(fraction+1.0));
    fprintf(stderr,"Mass 1 and 2 are %g and %g\n",m1,m2); 

    double eps = 1.0/(sqrt((double)1024));
    double m0 = 1.0/(double)1024;
    double eps1 = eps*cbrt(m1/m0);
    double eps2 = eps*cbrt(m2/m0);


    // Generate a particles distribution.
    mycount = 0;

    for(int i=0;i<Number;i++){
        double Pos[3];
        if(i%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeDM;
            Pbody[mycount]->GlobalID = i;

            if(i%2 == 0){
                Pbody[mycount]->Mass = m1;
                Pbody[mycount]->Eps = eps1;
            }else{
                Pbody[mycount]->Mass = m2;
                Pbody[mycount]->Eps = eps2;
            }

            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
            Pbody[mycount]->Pos[0] = Radius*Pos[0];
            Pbody[mycount]->Pos[1] = Radius*Pos[1];
            Pbody[mycount]->Pos[2] = Radius*Pos[2];

            Pbody[mycount]->Vel[0] = 
            Pbody[mycount]->Vel[1] = 
            Pbody[mycount]->Vel[2] = 0.0;

            Pbody[mycount]->InteractionList = 1;

            mycount ++;
        } else {
            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
        }
    }
    dprintlmpi(mycount);

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 1.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 10*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/CC.ASCII");
    strcpy(Pall.BaseFileName,"./data/CC");
    strcpy(Pall.RestartFileName,"./data/CC.dump");

    return;
}

void RestoreSphericalParticleDistributionErrorEstimation(const int Number, double Pos[][3], double Mass[], double Eps[]){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);
    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1;
    Pall.GravConst = 1;

    // Allocate Particle Data Structures.
    int mycount = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            mycount ++;

    Pall.Ntotal_t = Number;
    Pall.Ntotal = mycount;

    int AllocationSize = 1.2*mycount; 
    GenerateStructPbody(AllocationSize);
    // Allocate Particle Data Structures.
    
    // Generate a particles distribution.
    mycount = 0;
    for(int i=0;i<Number;i++){
        Pbody[mycount]->Active = ON;
        Pbody[mycount]->Use = ON;
        Pbody[mycount]->Type = TypeDM;
        Pbody[mycount]->GlobalID = i;

        Pbody[mycount]->Pos[0] = Pos[i][0];
        Pbody[mycount]->Pos[1] = Pos[i][1];
        Pbody[mycount]->Pos[2] = Pos[i][2];
        Pbody[mycount]->Mass = Mass[i];
        Pbody[mycount]->Eps = Eps[i];

        Pbody[mycount]->Vel[0] = 
        Pbody[mycount]->Vel[1] = 
        Pbody[mycount]->Vel[2] = 0.0;

        Pbody[mycount]->InteractionList = 1;

        mycount ++;
    }
    dprintlmpi(mycount);

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 1.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 10*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/CC.ASCII");
    strcpy(Pall.BaseFileName,"./data/CC");
    strcpy(Pall.RestartFileName,"./data/CC.dump");

    return;
}

static void CalcPotentialDirectForColdCollaspe(void){

    int NumProcs = MPIGetNumProcs();
    double r,r2,xij[3],EPS;
    MPI_Status  mpi_status;

    /* Self-Domain Potential */
    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Pot = 0.e0;
        for(int j=0;j<Pall.Ntotal;j++){
            xij[0] = Pbody[i]->Pos[0]-Pbody[j]->Pos[0];
            xij[1] = Pbody[i]->Pos[1]-Pbody[j]->Pos[1];
            xij[2] = Pbody[i]->Pos[2]-Pbody[j]->Pos[2];
            r2 = NORM2(xij);
            EPS = 0.5*(Pbody[i]->Eps+Pbody[j]->Eps);
            //if(r2 < SQ(TINY*EPS)) continue;

            r = sqrt(r2+SQ(EPS));
            Pbody[i]->Pot -= Pbody[j]->Mass/r;
        }
    }

    struct StructPotential{
        double Pos[3];
        double Eps;
        double Mass;
    } *PotentialExportSend,*PotentialExportRecv;

    /* Get Other Domains and Calculate Gravity */
    for(int i=0;i<NumProcs-1;i++){
        int SendThisTime,RecvThisTime;
        int SendAlready = 0;
        int RecvAlready = 0;
        int SendData = CommunicationTable[i].SendSize;
        int RecvData = CommunicationTable[i].RecvSize;
        while((SendData != 0)&&(RecvData != 0)){
            SendThisTime = CommunicationTable[i].SendSize-SendAlready;
            CheckSizeofBufferExportSendIndex(SendThisTime,sizeof(struct StructPotential),0);
            PotentialExportSend = BufferExportSend[0];

            int NActives = 0;
            while(NActives != SendThisTime){
                PotentialExportSend[NActives].Pos[0] = Pbody[NActives+SendAlready]->Pos[0];
                PotentialExportSend[NActives].Pos[1] = Pbody[NActives+SendAlready]->Pos[1];
                PotentialExportSend[NActives].Pos[2] = Pbody[NActives+SendAlready]->Pos[2];
                PotentialExportSend[NActives].Mass = Pbody[NActives+SendAlready]->Mass;
                PotentialExportSend[NActives].Eps = Pbody[NActives+SendAlready]->Eps;
                NActives ++;
            }
            RecvThisTime = CommunicationTable[i].RecvSize-RecvAlready;
            CheckSizeofBufferExportRecv(RecvThisTime,sizeof(struct StructPotential));
            PotentialExportRecv = BufferExportRecv;

            // Communication
            MPI_Sendrecv(PotentialExportSend,SendThisTime*sizeof(struct StructPotential),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_POTENTIAL_DIRECT_SENDRECV,
                PotentialExportRecv,RecvThisTime*sizeof(struct StructPotential),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_POTENTIAL_DIRECT_SENDRECV,
                        MPI_COMM_WORLD,&mpi_status);

            // Calculation
            for(int j=0;j<Pall.Ntotal;j++){
                for(int k=0;k<RecvThisTime;k++){
                    xij[0] = Pbody[j]->Pos[0]-PotentialExportRecv[k].Pos[0];
                    xij[1] = Pbody[j]->Pos[1]-PotentialExportRecv[k].Pos[1];
                    xij[2] = Pbody[j]->Pos[2]-PotentialExportRecv[k].Pos[2];
                    r2 = NORM2(xij);
                    EPS = 0.5*(Pbody[j]->Eps+PotentialExportRecv[k].Eps);
                    //if(r2 < SQ(TINY*EPS)) continue;

                    r = sqrt(r2+SQ(EPS));
                    Pbody[j]->Pot -= PotentialExportRecv[k].Mass/r;
                }
            }

            SendData -= SendThisTime;
            RecvData -= RecvThisTime;
            SendAlready += SendThisTime;
            RecvAlready += RecvThisTime;
        }
    }


    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Pot += Pbody[i]->Mass/Pbody[i]->Eps;
        Pbody[i]->Pot *= 0.5*Pall.GravConst*Pbody[i]->Mass;
    }

    if(NumProcs > 1){
        Delete(PotentialExportSend);
        Delete(PotentialExportRecv);
    }

    return;
}

/*
 * This function makes an initial condition for a cold collapse test with two
 * different mass particles. A half of particles have the mass of m1, and others
 * have the mass of m2. Since the total mass of the system is 1, the particle
 * masses should satisfy the relation: 1 = N(m1+m2)/2. When we adopt m1:m2 =
 * 1:2, the relation is 1 = 3 N m1/2. Thus m1 = 2/(3N). The third parameter of
 * this function gives m1/m2.
 */
void InitColdCollapseTestMixed(const int Number, const double rv, const double fraction, const int MakePotential){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);
    Pall.UnitTime = Pall.UnitLength = Pall.UnitMass = 1;
    //Pall.GravConst = GetUnitGravitationalConstant();
    Pall.GravConst = 1;

    // Allocate Particle Data Structures.
    int mycount = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            mycount ++;

    Pall.Ntotal_t = Number;
    Pall.Ntotal = mycount;

    int AllocationSize = mycount; 
    GenerateStructPbody(AllocationSize);
    // Allocate Particle Data Structures.
    
    
    double m1 = 2.0*fraction/((double)Number*(fraction+1.0));
    double m2 = 2.0/((double)Number*(fraction+1.0));
    fprintf(stderr,"Mass 1 and 2 are %g and %g\n",m1,m2); 

    double eps = 1.0/(sqrt((double)Number));
    double m0 = 1.0/(double)Number;
    double eps1 = eps*cbrt(m1/m0);
    double eps2 = eps*cbrt(m2/m0);

    // Generate a particles distribution.
    mycount = 0;


    for(int i=0;i<Number;i++){
        double Pos[3];
        if(i%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeDM;
            Pbody[mycount]->GlobalID = i;

            if(i%2 == 0){
                Pbody[mycount]->Mass = m1;
                Pbody[mycount]->Eps = eps1;
            }else{
                Pbody[mycount]->Mass = m2;
                Pbody[mycount]->Eps = eps2;
            }

            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
            Pbody[mycount]->Pos[0] = Pos[0];
            Pbody[mycount]->Pos[1] = Pos[1];
            Pbody[mycount]->Pos[2] = Pos[2];

            Pbody[mycount]->InteractionList = 1;

            mycount ++;
        } else {
            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
        }
    }
    dprintlmpi(mycount);

    double sigma;
    if(MakePotential == ON){
        //PotentialDirect();
        CalcPotentialDirectForColdCollaspe();

        // Calc Potential.
        double W = 0.e0;
        for(int i=0;i<mycount;i++)
            W += Pbody[i]->Pot;
        eprintlmpi(W);
        double GW;
        MPI_Allreduce(&W,&GW,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        W = GW;

        double Mtotal = 1.e0;
        sigma = sqrt(2.0*fabs(W)*rv/(3.0*Mtotal));

        eprintlmpi(sigma);
    } else {
        sigma = 0.e0;
    }

    // Calc sigma and add velocity for particles.
    mycount = 0;
    double Vel[3];
    for(int i=0;i<Pall.Ntotal;i++){
        if(i%NProcs == MyID){
            Pbody[mycount]->Vel[0] = sigma*Gaussian();
            Pbody[mycount]->Vel[1] = sigma*Gaussian();
            Pbody[mycount]->Vel[2] = sigma*Gaussian();
            mycount ++;
        } else {
            Vel[0] = gsl_ran_gaussian(RandomGenerator,sigma);
            Vel[1] = gsl_ran_gaussian(RandomGenerator,sigma);
            Vel[2] = gsl_ran_gaussian(RandomGenerator,sigma);
        }
    }

    double total_mass = 0.e0;
    double MeanVel[3] = {0.e0,0.e0,0.e0};
    for(int i=0;i<Pall.Ntotal;i++){
        MeanVel[0] += Pbody[i]->Mass*Pbody[i]->Vel[0];
        MeanVel[1] += Pbody[i]->Mass*Pbody[i]->Vel[1];
        MeanVel[2] += Pbody[i]->Mass*Pbody[i]->Vel[2];
        total_mass += Pbody[i]->Mass;
    }
    //gprintlmpi(total_mass);
    //fflush(NULL);
    MeanVel[0] /= total_mass;
    MeanVel[1] /= total_mass;
    MeanVel[2] /= total_mass;
    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Vel[0] -= MeanVel[0];
        Pbody[i]->Vel[1] -= MeanVel[1];
        Pbody[i]->Vel[2] -= MeanVel[2];
    }

    /*
    for(int i=0;i<Number;i++){
        fprintf(stderr,"%d %g %g %g %g %g %g %g %g\n",i,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Mass,Pbody[i]->Eps);
    }
    exit(1);
    */

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 1.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 10*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/CC.ASCII");
    strcpy(Pall.BaseFileName,"./data/CC");
    strcpy(Pall.RestartFileName,"./data/CC.dump");

    return;
}
#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR

void InitSphereWithMeanParticleDistance(const int Number){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MyID);

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    dprintlmpi(AllocationSize);
    
    PbodySize = AllocationSize;
    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(AllocationSize*sizeof(StructPbodyptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<AllocationSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Pbody[i] = PbodyElements+i;

    // Generate a particles distribution.
    int mycount = 0; double Pos[3];
    //double eps = 1.0/(sqrt((double)Number));

    double TotalVolume = (4.0*PI/3.0);
    double SubVolume = TotalVolume/(double)Number;
    double eps = cbrt(3.0*SubVolume/(4.0*PI));

    for(int i=0;i<Number;i++){
        if(i%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeDM;
            Pbody[mycount]->GlobalID = i;

            Pbody[mycount]->Mass = 1.e0/((double)Number);
            Pbody[mycount]->Eps = eps;

            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
            Pbody[mycount]->Pos[0] = Pos[0];
            Pbody[mycount]->Pos[1] = Pos[1];
            Pbody[mycount]->Pos[2] = Pos[2];

            mycount ++;
        } else {
            do {
                Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
                Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            } while ( NORM2(Pos) > 1.0 );
        }
    }
    dprintlmpi(mycount);

    Pall.Ntotal = mycount;
    Pall.Ntotal_t = Number;
    Pall.GravConst = 1.e0;


    Pall.Ntotal = mycount;
    Pall.Ntotal_t = Number;

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 10.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    InitLogFiles();

    return;
}

// Return total energy.
double CalcTotalEnergyForColdCollapse(void){

    double EGlobal;
    double ETotal = 0.e0;

    // get ETotal.
    for(int i=0;i<Pall.Ntotal;i++){
        ETotal += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel)+Pbody[i]->Pot;
    }

    MPI_Allreduce(&ETotal,&EGlobal,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    return EGlobal;
}

// cat ./InitColdCollapse.00 | sort -n -k 7 | awk '{print $7, NR}' >& s
void OutPutColdCollapse(char *fname){

    FILE *fp;
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%g %g %g %g %g %g %g %g\n",
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                NORM(Pbody[i]->Pos),NORM(Pbody[i]->Vel));
    }
    fclose(fp);

    return;
}

// Navarro White test.
static void UpdateVel(const double lambda);
static double AddSpin(void);
static void GetPotential(void);
static void SetParticles(const int NGrid);

void InitializeNavarroWhiteTest(const int NGrid){

    // allocate 
    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MPIGetMyID());

    // Set unit
    Pall.UnitLength = MPC_CGS;
    Pall.UnitTime = 10.0*GIGAYEAR_CGS;
    Pall.UnitMass = 1.e+11*MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    eprintlmpi(Pall.GravConst);
    eprintlmpi(Pall.DegreeOfFreedom);
    eprintlmpi(Pall.HeliumWeight);
    eprintlmpi(Pall.MeanMolecularWeight);
    eprintlmpi(Pall.ConvertUtoT);
    eprintlmpi(Pall.ConvertTtoU);

    eprintlmpi(Pall.ConvertDensityToCGS);
    eprintlmpi(Pall.ConvertNumberDensityToCGS);

    // Set particles
	SetParticles(NGrid);

	//GetPotential();
    PotentialDirect();
	double lambda = AddSpin();

	while(lambda < 0.0999 || 0.101 < lambda){
		UpdateVel(lambda);
		lambda = AddSpin();
	}

    //char fname[MaxCharactersInLine];
    //sprintf(fname,"NavarroWhiteInit.%02d",MPIGetMyID());
    //OutPutNavarroWhite(fname);

    char fname_gas[MaxCharactersInLine],
            fname_star[MaxCharactersInLine],fname_dm[MaxCharactersInLine];
    sprintf(fname_gas,"InitNavarroWhite.Hydro.%02d.%02d",MPIGetMyID(),MPIGetNumProcs());
    sprintf(fname_star,"InitNavarroWhite.Star.%02d.%02d",MPIGetMyID(),MPIGetNumProcs());
    sprintf(fname_dm,"InitNavarroWhite.DM.%02d.%02d",MPIGetMyID(),MPIGetNumProcs());
    OutPutNavarroWhite(fname_gas,fname_star,fname_dm);

    InitLogFiles();

	return ;
}

static void SetParticles(const int NGrid){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
	double r;

    double NX,NY,NZ; 
    NX = NY = NZ = NGrid;

    // count 
    double Pos[3];
	double dx = 2.e0/(double)(NX-1);
	int count = 0;
    int mycount = 0;
	for(int l=0;l<2;l++){
		for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
				for(int k=0;k<NZ;k++){
					Pos[0] = -1.e0+i*dx;
					Pos[1] = -1.e0+j*dx;
					Pos[2] = -1.e0+k*dx;

                    r = NORM(Pos);
                    if(r>TINY){
					    /*** make unit vector ***/
                        Pos[0] /= r;
                        Pos[1] /= r;
                        Pos[2] /= r;
                            
                        /*** strech grid(r_new = r_old**3) ***/
                        Pos[0] *= 0.1*(r)*sqrt(r);
                        Pos[1] *= 0.1*(r)*sqrt(r);
                        Pos[2] *= 0.1*(r)*sqrt(r);
				    } else {
                        Pos[0] = 0.e0;
                        Pos[1] = 0.e0;
                        Pos[2] = 0.e0;
                    }

					if(0.01>=NORM2(Pos)){
                        if(count%NProcs == MyID){
                            mycount ++;
                        }
						count ++;
					}
				}
			}
		}
	}

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = mycount;
    Pall.Ntotal_t = count;

    Pall.Nhydro = mycount/2;
    Pall.Nhydro_t = count/2;

    Pall.NDM = mycount/2;
    Pall.NDM_t = count/2;

    int AllocationSize = mycount; 
    int HydroAllocationSize = mycount/2; 
    PbodySize = mycount;
    PhydroSize = mycount/2;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(HydroAllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,HydroAllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);

    for(int i=0;i<HydroAllocationSize-1;i++)
        PhydroElements[i].Next = &(PhydroElements[i+1]);

    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[HydroAllocationSize-1].Next = NULL;


    for(int i=0;i<AllocationSize;i++)
        Pbody[i] = PbodyElements+i;

    for(int i=0;i<HydroAllocationSize;i++){
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


	dx = 2.e0/(double)(NX-1);
	count = mycount = 0;
	for(int l=0;l<2;l++){
		for(int i=0;i<NX;i++){
			for(int j=0;j<NY;j++){
				for(int k=0;k<NZ;k++){
					Pos[0] = -1.e0+i*dx;
					Pos[1] = -1.e0+j*dx;
					Pos[2] = -1.e0+k*dx;

					r = NORM(Pos);
                    if(r>TINY){
                        /*** make unit vector ***/
                        Pos[0] /= r;
                        Pos[1] /= r;
                        Pos[2] /= r;
				
                        /*** strech grid(r_new = r_old**3) ***/
                        Pos[0] *= 0.1*(r)*sqrt(r);
                        Pos[1] *= 0.1*(r)*sqrt(r);
                        Pos[2] *= 0.1*(r)*sqrt(r);
                    } else {
                        Pos[0] = 0.e0;
                        Pos[1] = 0.e0;
                        Pos[2] = 0.e0;
                    }
				
					if(0.01>=NORM2(Pos)){
                        if(count%NProcs == MyID){
                            Pbody[mycount]->Pos[0] = Pos[0];
                            Pbody[mycount]->Pos[1] = Pos[1];
                            Pbody[mycount]->Pos[2] = Pos[2];

                            Pbody[mycount]->Vel[0] = -Pos[1];
                            Pbody[mycount]->Vel[1] = +Pos[0];
                            Pbody[mycount]->Vel[2] = 0.e0;

                            Pbody[mycount]->GlobalID = count;
                            mycount ++;
                        }
						count ++;
					}
				}
			}
		}
	}
//#define E_CDM(x)    ( pow((x)/(1.0e-5)*CUBE(5.0e-4),1.0/3.0) )
//#define E_SPH(x)    ( pow((x)/(1.0e-5)*CUBE(5.0e-4),1.0/3.0) )

	//double imass = 0.1/((double)count/2); 
	double imass = 10.0/((double)count/2); 
    double EpsHydro = cbrt((0.1*imass)/(1.0e-5)*CUBE(5.0e-4));
    double EpsDM    = EpsHydro;
    //double EpsDM    = cbrt((0.9*imass)/(1.0e-5)*CUBE(5.0e-4));

//#define ConversionFactorTtoU    (2.1619e-6)
    ///dprintlmpi(mycount);
    ///dprintlmpi(count);
    ///dprintlmpi(mycount/2);
    ///dprintlmpi(Pall.Ntotal);
    ///dprintlmpi(Pall.NDM);
    ///dprintlmpi(Pall.Nhydro);

    double TotalMass = 0.0;

    Pall.Ntotal = 0;
    Pall.Nhydro = 0;
    Pall.NDM = 0;
	int hcount = mycount/2;
	for(int l=0;l<2;l++){
		for(int i=0;i<hcount;i++){
			if(l == 0){
                Pbody[i]->Active = ON;
                Pbody[i]->Use = ON;
                PbodyHydro(i)->Use = ON;
                Pbody[i]->Type = TypeHydro;

                Pbody[i]->Mass = 0.1*imass;
#if (UseSFModelSpawn) 
                PbodyHydro(i)->SpawnMass = Pbody[i]->Mass/((double)MaxSpawnTimes);
#endif
                Pbody[i]->Eps = EpsHydro;
                PbodyHydroKernel(i) = EpsHydro;
                PbodyHydroU(i) = 1000.e0*Pall.ConvertTtoU;
                TotalMass += 0.1*imass;
                Pall.Nhydro ++;
			}else{
                Pbody[i+hcount]->Active = ON;
                Pbody[i+hcount]->Use = ON;
                Pbody[i+hcount]->Type = TypeDM;

                Pbody[i+hcount]->Mass = 0.9*imass;
                Pbody[i+hcount]->Eps = EpsDM;
                TotalMass += 0.9*imass;
                Pall.NDM ++;
			}
            Pall.Ntotal ++;
		}
	}
    dlprintlmpi(Pall.Ntotal);
    dlprintlmpi(Pall.NDM);
    dlprintlmpi(Pall.Nhydro);
    eprintlmpi(TotalMass);

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }


    Pall.RunStatus = NewSimulation;

    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = 0.8;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
            Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.AdaptiveSofteningFactor = 1.e0;

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);


#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 0.75;
    Pall.ViscousS = 0.75;
    Pall.ViscousL = 0.75;
#endif // USE_VARIABLE_ALPHA

    strcpy(Pall.ASCIIFileName,"./data/NFW.ASCII");
    strcpy(Pall.BaseFileName,"./data/NFW");
    strcpy(Pall.RestartFileName,"./data/NFW.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.01;

    FileOutPutConstantInterval();
    InitializeRandomGenerator(1977);

	return ;

}

// Need Parallel Version
static void GetPotential(void){

	int i,k;
	double r,r2;

//#pragma omp parallel for private(k,r)  schedule(guided,1)
	for(i=0;i<Pall.Ntotal;i++){
		Pbody[i]->Pot = 0.e0;
		for(k=0;k<Pall.Ntotal;k++){
			r2 = DISTANCE2(Pbody[i]->Pos,Pbody[k]->Pos);
			if( (r2 > SQ(0.5*(Pbody[i]->Eps+Pbody[k]->Eps)*TINY)) ){
				r2 += SQ(0.5*(Pbody[i]->Eps+Pbody[k]->Eps));
				r = sqrt(r2);
				Pbody[i]->Pot += Pbody[k]->Mass/r;
			}
		}
		Pbody[i]->Pot *= -0.5*Pall.GravConst*Pbody[i]->Mass;
	}

	return;
}

static double AddSpin(void){

	int i;
	double Etotal,Ek,Ep,Mtotal,xv[3];

	Etotal = Mtotal = Ek = Ep = xv[0] = xv[1] = xv[2] = 0.e0;
	for(i=0;i<Pall.Ntotal;i++){
		Ek += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel);
        Ep += Pbody[i]->Pot;
		Etotal += 0.5*Pbody[i]->Mass*NORM2(Pbody[i]->Vel) + Pbody[i]->Pot;

		xv[0] += Pbody[i]->Mass*(Pbody[i]->Pos[1]*Pbody[i]->Vel[2] - Pbody[i]->Pos[2]*Pbody[i]->Vel[1]);
		xv[1] += Pbody[i]->Mass*(Pbody[i]->Pos[2]*Pbody[i]->Vel[0] - Pbody[i]->Pos[0]*Pbody[i]->Vel[2]);
		xv[2] += Pbody[i]->Mass*(Pbody[i]->Pos[0]*Pbody[i]->Vel[1] - Pbody[i]->Pos[1]*Pbody[i]->Vel[0]);

		Mtotal += Pbody[i]->Mass;
	}

    // communication!
    double AMEnergyMass[5],GlobalAMEnergyMass[5];
    AMEnergyMass[0] = xv[0];
    AMEnergyMass[1] = xv[1];
    AMEnergyMass[2] = xv[2];
    AMEnergyMass[3] = Etotal;
    AMEnergyMass[4] = Mtotal;
    MPI_Allreduce(AMEnergyMass,GlobalAMEnergyMass,5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    xv[0] = GlobalAMEnergyMass[0];
    xv[1] = GlobalAMEnergyMass[1];
    xv[2] = GlobalAMEnergyMass[2];
    Etotal = GlobalAMEnergyMass[3];
    Mtotal = GlobalAMEnergyMass[4];
    // end communication!

	double AMomentum = NORM(xv); // not specific.

	fprintf(stderr,"(Etotal,AMomentum,Mass) = (%g ,%g ,%g)\n",Etotal,AMomentum,Mtotal);

	double lambda = AMomentum*sqrt(fabs(Etotal))/(Pall.GravConst*SQ(Mtotal)*sqrt(Mtotal));
	fprintf(stderr,"Lambda = %g : in function %s\n",lambda,__FUNCTION__);

	return ( lambda );
}

static void UpdateVel(const double lambda){

	double factor = sqrt(sqrt(0.1/lambda));
	for(int i=0;i<Pall.Ntotal;i++){
		Pbody[i]->Vel[0] *= factor;
		Pbody[i]->Vel[1] *= factor;
		Pbody[i]->Vel[2] *= factor;
	}

	return ;
}

void OutPutDarkMatter(char *fname_dm){

	FILE *fp_dm;

    FileOpen(fp_dm,fname_dm,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            fprintf(fp_dm,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                    Pbody[i]->Mass);
        }
    }
    fclose(fp_dm);

	return;
}

void OutPutNavarroWhite(char *fname_hydro, char *fname_star, char *fname_dm){

	FILE *fp_hydro,*fp_star,*fp_dm;

    FileOpen(fp_hydro,fname_hydro,"w");
    FileOpen(fp_star,fname_star,"w");
    FileOpen(fp_dm,fname_dm,"w");

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
            fprintf(fp_hydro,"%ld %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                    Pbody[i]->Mass,Pall.ConvertNumberDensityToCGS*PbodyHydro(i)->Rho,
                    PbodyHydroU(i)*Pall.ConvertUtoT,PbodyHydroKernel(i));
        } else if(Pbody[i]->Type == TypeStar){
            fprintf(fp_star,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                    Pbody[i]->Mass);
        } else {
            fprintf(fp_dm,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                    Pbody[i]->Mass);
        }
    }

    fclose(fp_hydro);
    fclose(fp_star);
    fclose(fp_dm);

	return;
}


#if 0
static void CircularProfile(void){

	int i;
	FILE *fp;
	
	FileOpen(fp,"circular.data","w");
	for(i=0;i<count;i++)
		fprintf(fp,"%e %e %e %e %e %e %e\n",
			NORM(x[i]),sqrt(SQ(x[i][0])+SQ(x[i][1])),NORM(v[i]),sqrt(pot[i]/m[i]),
			sqrt(G*10.e0*NORM2(x[i])/CUBE(0.1)),10.e0*CUBE(NORM(x[i]))/CUBE(0.1),
			NORM(v[i])/sqrt(SQ(x[i][0])+SQ(x[i][1])));
	fclose(fp);

	return;
}
#endif

void InitTreeTest(const int Number){
    
    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MyID);

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Number;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    dprintlmpi(AllocationSize);
    
    PbodySize = AllocationSize;
    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(AllocationSize*sizeof(StructPbodyptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<AllocationSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Pbody[i] = PbodyElements+i;
    // Allocate Particle Data Structures.

    // Generate a particles distribution.
    int mycount = 0; double Pos[3];
    double eps = 1.0/(sqrt((double)Number));
    for(int i=0;i<Number;i++){
        if(i%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeDM;
            Pbody[mycount]->GlobalID = i;

            Pbody[mycount]->Mass = 1.e0/((double)Number);
            Pbody[mycount]->Eps = eps;

            Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;

            Pbody[mycount]->Pos[0] = Pos[0];
            Pbody[mycount]->Pos[1] = Pos[1];
            Pbody[mycount]->Pos[2] = Pos[2];

            mycount ++;
        } else {
            Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        }
    }

    Pall.Ntotal = mycount;
    Pall.Ntotal_t = Number;
    Pall.NActives = mycount;
    Pall.NActives_t = Number;
    Pall.GravConst = 1.e0;
    //PotentialDirect();

    // Calc Potential.
    double W = 0.e0;
    for(int i=0;i<mycount;i++)
        W += Pbody[i]->Pot;
    eprintlmpi(W);
    double GW;
    MPI_Allreduce(&W,&GW,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    W = GW;

    //double Mtotal = 1.e0;
    //double sigma = sqrt(2.0*fabs(W)*rv/(3.0*Mtotal));
    double sigma = 1.e0;

    // Calc sigma and add velocity for particles.
    mycount = 0;
    double Vel[3];
    for(int i=0;i<Number;i++){
        if(i%NProcs == MyID){
            Pbody[mycount]->Vel[0] = sigma * Gaussian();
            Pbody[mycount]->Vel[1] = sigma * Gaussian();
            Pbody[mycount]->Vel[2] = sigma * Gaussian();
            mycount ++;
        } else {
            Vel[0] = sigma * Gaussian();
            Vel[1] = sigma * Gaussian();
            Vel[2] = sigma * Gaussian();
        }
    }

    Pall.Ntotal = mycount;
    Pall.Ntotal_t = Number;

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 1.e0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    InitLogFiles();

    return;
}


// Isothermal Spherical Collapse test.
void InitIsothermalSphericalCollapseTest(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = number;

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));

    // Set unit
    Pall.UnitLength = 7.0e+16;
    Pall.UnitTime = 1.7888e+12; // sqrt(3*M_PI/(32*G*rho_c))
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 5.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 2.0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    eprintlmpi(Pall.GravConst);
    eprintlmpi(Pall.DegreeOfFreedom);
    eprintlmpi(Pall.HeliumWeight);
    eprintlmpi(Pall.MeanMolecularWeight);
    eprintlmpi(Pall.ConvertUtoT);
    eprintlmpi(Pall.ConvertTtoU);

    eprintlmpi(Pall.ConvertDensityToCGS);
    eprintlmpi(Pall.ConvertNumberDensityToCGS);

    eprintlmpi(GetUnitVel());

    // Set Parameters
    // Rsphere, MCloud, cs, and dOmega/dt.
    double Rsphere = 7.e+16/Pall.UnitLength;
    double MCloud = 1.e0; // 1 Msun
    double Temperature = 10.0; // 10K
    double cs = sqrt(UNIVERSAL_GAS_CONSTANT_CGS*AVOGADROS_NUMBER_CGS*Temperature/
            (Pall.MeanMolecularWeight*PROTON_MASS_CGS))*GetUnitVel(); // 0.2 km/s ok?
    //double cs = sqrt(Pall.DegreeOfFreedom*UNIVERSAL_GAS_CONSTANT_CGS*AVOGADROS_NUMBER_CGS*Temperature/
            //(Pall.MeanMolecularWeight*PROTON_MASS_CGS))*GetUnitVel(); // 0.32 km/s ok?
    double omega = 3.04e-13*Pall.UnitTime;

    eprintlmpi(cs);

    double Nhalf_real = (double)NX;

    double Pos[3];
    int count = 0;
    int mycount = 0;
    for(int i=-NX;i<NX;i++){
    for(int j=-NY;j<NY;j++){
    for(int k=-NZ;k<NZ;k++){
        Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
        Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
        Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

        if(Rsphere>NORM(Pos)){
            if(count%NProcs == MyID){
                mycount ++;
            }
            count ++;
        }
    }
    }
    }


    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount;
    PbodySize = mycount;
    PhydroSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = MCloud/(double)count;
    double eps = 1.e-5*PC_CGS/Pall.UnitLength;
    //double Uinit = SQ(cs);

    count = 0; mycount = 0;
    for(int i=-NX;i<NX;i++){
    for(int j=-NY;j<NY;j++){
    for(int k=-NZ;k<NZ;k++){
        Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
        Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
        Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

        if(Rsphere>NORM(Pos)){
            if(count%NProcs == MyID){
                Pbody[mycount]->Active = ON;
                Pbody[mycount]->Use = ON;
                Pbody[mycount]->Type = TypeHydro;
                Pbody[mycount]->GlobalID = count;

                Pbody[mycount]->Pos[0] = Pos[0];
                Pbody[mycount]->Pos[1] = Pos[1];
                Pbody[mycount]->Pos[2] = Pos[2];

                double Rxy = sqrt(SQ(Pos[0])+SQ(Pos[1]));

                if(Rxy < TINY*Rsphere)
                    exit(444);

                Pbody[mycount]->Mass = mass;

                //v_tangential = (Rxy * KM_PER_PC) * omega_in_rad_per_sec;
                double v_tangential = Rxy * omega;
                double vx=  -v_tangential * Pos[1]/Rxy;
                double vy=   v_tangential * Pos[0]/Rxy;

                Pbody[mycount]->Vel[0] = vx;
                Pbody[mycount]->Vel[1] = vy;
                Pbody[mycount]->Vel[2] = 0;

                Pbody[mycount]->Eps = eps;

                PbodyHydro(mycount)->Use = ON;
                PbodyHydro(mycount)->Kernel = 2.0*Rsphere/NX;

                //PbodyHydro(mycount)->U = Uinit;
                PbodyHydro(mycount)->U = 0.e0;

                mycount ++;
            }
            count ++;
        }
    }
    }
    }

    double Tpot = Pall.GravConst*MCloud/Rsphere;
    double Trot = (3.0/10.0)*MCloud*SQ(Rsphere)*SQ(omega); 
    double Tth = MCloud*GetUnitSpecificEnergy()*
        (Pall.DegreeOfFreedom*BOLTZMANN_CONSTANT_CGS*Temperature/
         (Pall.MeanMolecularWeight*PROTON_MASS_CGS));

    fprintf(stderr,"Tpot = %g, Trot = %g, Tth = %g, alpha = %g(0.52), beta = %g(0.08)\n",
            Tpot,Trot,Tth,Tth/Tpot,Trot/Tpot);

    /*
    double Mt = 0.e0; double Rt = 0.e0; double Tt = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        Mt += PhydroMass(i);
        Rt += PhydroMass(i)*NORM2(PhydroPos(i))*NORM(PhydroVel(i));
        Tt += PhydroMass(i)*Phydro[i]->U;
    }
    */
    // Check Total Mass

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"InitIsothermal.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    exit(0);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.CS = cs;

    Pall.TEnd = 1.2;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
            Pall.TEnd,Pall.TEnd*Pall.UnitTime);

    Pall.Gamma = 1.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.AdaptiveSofteningFactor = 1.e0;

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    strcpy(Pall.ASCIIFileName,"./data/Isothermal.ASCII");
    strcpy(Pall.BaseFileName,"./data/Isothermal");
    strcpy(Pall.RestartFileName,"./data/Isothermal.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1;

    OutPutASCIIDATA();

    FileOutPutConstantInterval();
    InitializeRandomGenerator(1977);

    InitLogFiles();


    return;
}

/*
 * This function generates the initial condition for the standard isothermal
 * test run done by e.g., Boss and Bodenheimer (1979) and Burkert and
 * Bodenheimer (1993).
 */
void InitStandardIsothermalTestCase(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = number;

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    // Set unit
    Pall.UnitLength = 3.2e+16;
    Pall.UnitTime = 5.52e+11; // sqrt(3*M_PI/(32*G*rho_c))
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 5.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 2.0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    eprintlmpi(Pall.GravConst);
    eprintlmpi(Pall.DegreeOfFreedom);
    eprintlmpi(Pall.HeliumWeight);
    eprintlmpi(Pall.MeanMolecularWeight);
    eprintlmpi(Pall.ConvertUtoT);
    eprintlmpi(Pall.ConvertTtoU);

    eprintlmpi(Pall.ConvertDensityToCGS);
    eprintlmpi(Pall.ConvertNumberDensityToCGS);

    eprintlmpi(GetUnitVel());

    // Set Parameters
    // Rsphere, MCloud, cs, and dOmega/dt.
    double Rsphere = 1.e0;
    double MCloud = 1.e0; // 1 Msun
    double Temperature = 10.0; // 10K
    double cs = sqrt(UNIVERSAL_GAS_CONSTANT_CGS*AVOGADROS_NUMBER_CGS*Temperature/
            (Pall.MeanMolecularWeight*PROTON_MASS_CGS))*GetUnitVel(); // 0.2 km/s ok?
    //double cs = sqrt(Pall.DegreeOfFreedom*UNIVERSAL_GAS_CONSTANT_CGS*AVOGADROS_NUMBER_CGS*Temperature/
            //(Pall.MeanMolecularWeight*PROTON_MASS_CGS))*GetUnitVel(); // 0.32 km/s ok?
    double omega = 1.6e-12*Pall.UnitTime;

    eprintlmpi(cs);
    fprintf(stderr,"Sound Speed = %g [cm/s]\n",
            sqrt(UNIVERSAL_GAS_CONSTANT_CGS*AVOGADROS_NUMBER_CGS*Temperature/
                (Pall.MeanMolecularWeight*PROTON_MASS_CGS)));

    double Nhalf_real = (double)NX;

    double Pos[3];
    int count = 0;
    int mycount = 0;
    for(int i=-NX;i<NX;i++){
    for(int j=-NY;j<NY;j++){
    for(int k=-NZ;k<NZ;k++){
        Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
        Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
        Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

        if(Rsphere>NORM(Pos)){
            if(count%NProcs == MyID){
                mycount ++;
            }
            count ++;
        }
    }
    }
    }


    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount;
    PbodySize = mycount;
    PhydroSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = MCloud/(double)count;
    double eps = 1.e-5*PC_CGS/Pall.UnitLength;
    //double Uinit = SQ(cs);

    double rho_crit = 7.e-24*SQ(Pall.Nhydro_t);
    eps = cbrt((6*50*mass/(M_PI*rho_crit*GetUnitDensity())))/4.0;
    fprintf(stderr,"esp = %g [pc], critical density = %g [g/cm^3]\n",eps,rho_crit);

    count = 0; mycount = 0;
    for(int i=-NX;i<NX;i++){
    for(int j=-NY;j<NY;j++){
    for(int k=-NZ;k<NZ;k++){
        Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
        Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
        Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

        if(Rsphere>NORM(Pos)){
            if(count%NProcs == MyID){
                Pbody[mycount]->Active = ON;
                Pbody[mycount]->Use = ON;
                Pbody[mycount]->Type = TypeHydro;
                Pbody[mycount]->GlobalID = count;

                Pbody[mycount]->Pos[0] = Pos[0];
                Pbody[mycount]->Pos[1] = Pos[1];
                Pbody[mycount]->Pos[2] = Pos[2];

                double Rxy = sqrt(SQ(Pos[0])+SQ(Pos[1]));

                if(Rxy < TINY*Rsphere)
                    exit(444);

                double phi;
                if(Pos[1] > 0.e0)
                    phi = acos(Pos[0]/Rxy);
                else
                    phi = 2.0*PI-acos(Pos[0]/Rxy);

                Pbody[mycount]->Mass = mass*(1.0+0.5*cos(2.0*phi));
                //Pbody[mycount]->Mass = mass;

                //v_tangential = (Rxy * KM_PER_PC) * omega_in_rad_per_sec;
                double v_tangential = Rxy * omega;
                double vx=  -v_tangential * Pos[1]/Rxy;
                double vy=   v_tangential * Pos[0]/Rxy;

                Pbody[mycount]->Vel[0] = vx;
                Pbody[mycount]->Vel[1] = vy;
                Pbody[mycount]->Vel[2] = 0;

                Pbody[mycount]->Eps = eps;

                PbodyHydro(mycount)->Use = ON;
                PbodyHydro(mycount)->Kernel = 2.0*Rsphere/NX;

                //PbodyHydro(mycount)->U = Uinit;
                PbodyHydro(mycount)->U = 0.e0;

                mycount ++;
            }
            count ++;
        }
    }
    }
    }

    double Tpot = Pall.GravConst*MCloud/Rsphere;
    double Trot = (3.0/10.0)*MCloud*SQ(Rsphere)*SQ(omega); 
    double Tth = MCloud*GetUnitSpecificEnergy()*
        (Pall.DegreeOfFreedom*BOLTZMANN_CONSTANT_CGS*Temperature/
         (Pall.MeanMolecularWeight*PROTON_MASS_CGS));

    fprintf(stderr,"Tpot = %g, Trot = %g, Tth = %g, alpha = %g(0.25), beta = %g(0.2)\n",
            Tpot,Trot,Tth,Tth/Tpot,Trot/Tpot);
    //exit(0);
    /*
    double Mt = 0.e0; double Rt = 0.e0; double Tt = 0.e0;
    for(int i=0;i<Pall.Nhydro;i++){
        Mt += PhydroMass(i);
        Rt += PhydroMass(i)*NORM2(PhydroPos(i))*NORM(PhydroVel(i));
        Tt += PhydroMass(i)*Phydro[i]->U;
    }
    */
    // Check Total Mass

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"InitIsothermal.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    exit(0);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.CS = cs;

    Pall.TEnd = 1.3;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
            Pall.TEnd,Pall.TEnd*Pall.UnitTime);

    Pall.Gamma = 1.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.AdaptiveSofteningFactor = 1.e0;

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    strcpy(Pall.ASCIIFileName,"./data/Isothermal.ASCII");
    strcpy(Pall.BaseFileName,"./data/Isothermal");
    strcpy(Pall.RestartFileName,"./data/Isothermal.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1;

    OutPutASCIIDATA();

    FileOutPutConstantInterval();

    InitLogFiles();

    return;
}

struct StructISC{
    double r;
    double sam;
    double mass;
} *ISC;

static int ISCCmp(const void *x, const void *y){
    const struct StructISC *pointer1 = x;
    const struct StructISC *pointer2 = y;
    if( pointer1->sam > pointer2->sam)
        return 1;
    else if( pointer1->sam < pointer2->sam)
        return -1;
    else
        return 0;
}

void FileOutputIsothermalSphericalCollapse(void){

    static double IOTiming[] = {0.0,0.8,1.1,1.2};
    int Num = sizeof(IOTiming)/sizeof(double);
    static int IOstep = 0;
    char fname[MaxCharactersInLine];
    FILE *fp;

    if(IOstep < Num){
        if(IOTiming[IOstep] < Pall.TCurrent){

            struct StructISC *ISC;
            ISC = malloc(sizeof(struct StructISC)*Pall.Nhydro);

            sprintf(fname,"IsothermalCollapse.%02d.%02d.%02d",MPIGetNumProcs(),MPIGetMyID(),IOstep);
            FileOpen(fp,fname,"w");
            for(int i=0;i<Pall.Nhydro;i++){
                ISC[i].r = NORM(PhydroPosP(i));
                ISC[i].sam = SQ(Pall.UnitLength)/Pall.UnitTime*(PhydroPosP(i)[0]*PhydroVel(i)[1]-PhydroPosP(i)[1]*PhydroVel(i)[0]);
                ISC[i].mass = Pall.UnitMass*PhydroMass(i);
            }

            qsort(ISC,Pall.Nhydro,sizeof(struct StructISC),(int(*)(const void*, const void*))ISCCmp);

            double Mass = 0.e0;
            for(int i=0;i<Pall.Nhydro;i++){
                Mass += ISC[i].mass;
               fprintf(fp,"%g %g %g\n",ISC[i].sam,Mass,ISC[i].r); 
            }
            fclose(fp);

            free(ISC);

            IOstep ++;
        }
    }

    return;
}


#ifdef TASK_M2_COLLAPSE
/* 
 * This function generates the initial condition of the isotheraml cloud
 * collapse test involving M2 perturbation. See Bate and Burkert (1997),
 * Truelove et al. (1997), and Matsumoto (2007).
 */
void InitM2SphericalCollapseTest(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = number;

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    // Set unit
    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = YEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    eprintlmpi(Pall.GravConst);
    eprintlmpi(Pall.DegreeOfFreedom);
    eprintlmpi(Pall.HeliumWeight);
    eprintlmpi(Pall.MeanMolecularWeight);
    eprintlmpi(Pall.ConvertUtoT);
    eprintlmpi(Pall.ConvertTtoU);

    eprintlmpi(Pall.ConvertDensityToCGS);
    eprintlmpi(Pall.ConvertNumberDensityToCGS);

    eprintlmpi(GetUnitVel());

    // Set Parameters
    // Rsphere, MCloud, cs, and dOmega/dt.
    double Rsphere = 5e+16/Pall.UnitLength; // 0.02pc
    double MCloud = 1.e0; // 1 Msun
    double cs = 0.1666*VELOCITY_KMS_CGS*GetUnitVel(); // 0.17 km/s
    //double omega = 7.2e-13*Pall.UnitTime; // 7.2x10^{-13} rad/s
    //double omega = 7.2e-13*GetUnitVel();
    double omega = 7.2e-13*Pall.UnitTime; // Rad/UnitTime

    eprintlmpi(cs);

    double Nhalf_real = (double)NX;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=-NX;i<NX;i++){
		for(int j=-NY;j<NY;j++){
			for(int k=-NZ;k<NZ;k++){
				Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
				Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
				Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

				if(Rsphere>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}


    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;

    int AllocationSize = mycount; 

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);

    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = MCloud/(double)count;
    //double eps = 0.1*(cbrt((double)1736/(double)count));
    //double eps = 1.e-5*PC_CGS/Pall.UnitLength;
    double eps = 10*AU_CGS/Pall.UnitLength;
    //double Uinit = SQ(cs)/(5.0/3.0 * (5.0/3.0-1.0));

    double rho_crit = 7.e-24*SQ(Pall.Nhydro_t);
    //eps = cbrt((6*50*mass/(M_PI*rho_crit*GetUnitDensity())))/5.6;
    eps = cbrt((6*50*mass/(M_PI*rho_crit*GetUnitDensity())))/4.0;
    fprintf(stderr,"esp = %g [pc], critical density = %g [g/cm^3]\n",eps,rho_crit);
    //exit(0);

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;

    double v_tangential,vx,vy;
    double Rxy;
    double phi;
    //double Uinit = 1.e-5*10.0*Pall.ConvertTtoU;
    double Uinit = SQ(cs)/Pall.GGm1;
	for(int i=-NX;i<NX;i++){
    for(int j=-NY;j<NY;j++){
    for(int k=-NZ;k<NZ;k++){
        Pos[0] = Rsphere*(i+0.5)/Nhalf_real;
        Pos[1] = Rsphere*(j+0.5)/Nhalf_real;
        Pos[2] = Rsphere*(k+0.5)/Nhalf_real;

        if(Rsphere>NORM(Pos)){
            if(count%NProcs == MyID){
                Pbody[mycount]->Active = ON;
                Pbody[mycount]->Use = ON;
                Pbody[mycount]->Type = TypeHydro;
                Pbody[mycount]->GlobalID = count;

                Pbody[mycount]->PosP[0] = Pbody[mycount]->Pos[0] = Pos[0];
                Pbody[mycount]->PosP[1] = Pbody[mycount]->Pos[1] = Pos[1];
                Pbody[mycount]->PosP[2] = Pbody[mycount]->Pos[2] = Pos[2];

                Rxy = sqrt(SQ(Pos[0])+SQ(Pos[1]));

                if(Rxy < TINY*Rsphere)
                    exit(444);

                if(Pos[1] > 0.e0)
                    phi = acos(Pos[0]/Rxy);
                else
                    phi = 2.0*PI-acos(Pos[0]/Rxy);

                Pbody[mycount]->Mass = mass*(1.0+0.1*cos(2.0*phi));

                //v_tangential = (Rxy * KM_PER_PC) * omega_in_rad_per_sec;
                v_tangential = Rxy * omega;
                vx=  -v_tangential * Pos[1]/Rxy;
                vy=   v_tangential * Pos[0]/Rxy;

                Pbody[mycount]->Vel[0] = vx;
                Pbody[mycount]->Vel[1] = vy;
                Pbody[mycount]->Vel[2] = 0;

                Pbody[mycount]->Eps = eps;

                PbodyHydro(mycount)->Use = ON;
                //PbodyHydro(mycount)->Kernel = 0.1*dx;
                PbodyHydro(mycount)->Kernel = Rsphere/NX;
                //eprintlmpi(PbodyHydro(mycount)->Kernel);

                //cout<<"Internal energy..."<<endl;
                // P = (gamma-1) rho u = cs_2 rho
                // thus 2/3 u = kT/mp = cs_2
                // thus  cs^2 = 3/2 kT/mp = u ; internal energy per mass

                //for(i=0;i<N_gas;i++)
                    //buf1[i]=pow(cs_in_cm_per_sec, 2) * (UnitMass/UnitEnergy);

                PbodyHydro(mycount)->U = Uinit;

#ifdef USE_CELIB //{
                //CELibSetPrimordialMetallicity(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements);
                CELibSetSolarMetallicity(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements);
                // CELibSetMetallicityWithSolarAbundancePattern(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements,
                // 0.5*CELibGetMetalFractionForSolarChemicalComposision());
                // CELibSetMetallicityWithSolarAbundancePattern(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements,0.02);
                double MassLightElements = Phydro[mycount]->Elements[CELibYield_H]+Phydro[mycount]->Elements[CELibYield_He];
                double Z = (Pbody[mycount]->Mass-MassLightElements)/Pbody[mycount]->Mass;
                PbodyHydro(mycount)->Z = PbodyHydro(mycount)->ZII = PbodyHydro(mycount)->ZIa = Z;
#else 
                PbodyHydro(mycount)->Z   = 0.02;
                PbodyHydro(mycount)->ZII = 0.02;
                PbodyHydro(mycount)->ZIa = 0.00;
#endif // USE_CELIB //}



                mycount ++;
            }
            count ++;
        }
    }
    }
	}

    ActivateAllparticles();

#if 0
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    ReConnectPointers();
    UpdateTotalNumber();
#endif

    // Check Total Mass

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"M2Collapse.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    fflush(NULL);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat M2Collapse.Init.* | sort -n > M2Collapse.dat");
        fflush(NULL);
        system("rm -rf ./M2Collapse.Init.*");
        fflush(NULL);
    }
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    // Pall.Ns = 64;
    // Pall.Npm = 2;
    Pall.Ns = 128;
    Pall.Npm = 8;
    Pall.CS = cs;

    //Pall.TEnd = 3.0; // ~1.3 dynamical time
    Pall.TEnd = 1.4*1.0774e+12/Pall.UnitTime;
    //Pall.TEnd = 3*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
            Pall.TEnd,Pall.TEnd*Pall.UnitTime);

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = 0.05*1.e+12/Pall.UnitTime;
    Pall.OutPutInterval = Pall.TEnd/100.0;
    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/M2.ASCII");
    strcpy(Pall.BaseFileName,"./data/M2");
    strcpy(Pall.RestartFileName,"./data/M2.dump");

    //FileOutPutConstantInterval();
    //InitLogFiles();
    //fflush(NULL);

    return;
}
#endif

#ifdef TASK_TURBULENCE
void InitTurbulentCloud(char fname[], const int mode, const double LeftMetal, const double RightMetal){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    // Set unit
    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    //double cs = 0.1666*VELOCITY_KMS_CGS*GetUnitVel(); // 0.17 km/s
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    // check the number of particles.

    FILE *fp;
    FileOpen(fp,fname,"rb");
    int collision,count;
    fread(&collision, sizeof(int), 1, fp);
    fread(&count, sizeof(int), 1, fp);
    fclose(fp);
    fprintf(stderr, "Collision flag = %d, Number of particles = %d\n",
            collision,count);

    int mycount = 0;
    for(int i=0;i<count;i++){
        if(i%NProcs == MyID){
            mycount ++;
        }
    }

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount;
    PbodySize = mycount;
    PhydroSize = mycount;

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);


    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    //double eps = 0.01*PC_CGS/Pall.UnitLength;
    // double eps = 5*AU_CGS/Pall.UnitLength;
    double eps = 50*AU_CGS/Pall.UnitLength;
    double kernel = 10*(PC_CGS/Pall.UnitLength)/cbrt(Pall.Ntotal_t);
    double Uinit = 10*Pall.ConvertTtoU;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"eps = %g, kernel = %g, Uinit = %g\n",
                eps,kernel,Uinit);

    FileOpen(fp,fname,"rb");

    for(int k=0;k<NProcs;k++){
        if(k == MyID){
            int dummy1;
            fread(&dummy1, sizeof(int), 1, fp);
            fread(&dummy1, sizeof(int), 1, fp);
            double dummy2;
            fread(&dummy2, sizeof(double), 1, fp);
            fread(&dummy2, sizeof(double), 1, fp);
            fread(&dummy2, sizeof(double), 1, fp);

            mycount = 0;
            for(int i=0;i<count;i++){
                struct StructPosVelMass{
                    double Pos[3];
                    double Vel[3];
                    double Mass;
                } PosVelMass;

                fread(&PosVelMass, sizeof(struct StructPosVelMass), 1, fp);

                if(i%NProcs == MyID){

                    Pbody[mycount]->Active = ON;
                    Pbody[mycount]->Use = ON;
                    Pbody[mycount]->Type = TypeHydro;
                    Pbody[mycount]->GlobalID = i;

                    Pbody[mycount]->Pos[0] = PosVelMass.Pos[0];
                    Pbody[mycount]->Pos[1] = PosVelMass.Pos[1];
                    Pbody[mycount]->Pos[2] = PosVelMass.Pos[2];

                    Pbody[mycount]->Mass = PosVelMass.Mass;

                    Pbody[mycount]->Vel[0] = PosVelMass.Vel[0];
                    Pbody[mycount]->Vel[1] = PosVelMass.Vel[1];
                    Pbody[mycount]->Vel[2] = PosVelMass.Vel[2];

                    Pbody[mycount]->Eps = eps;

                    PbodyHydro(mycount)->Use = ON;
                    PbodyHydro(mycount)->Kernel = kernel;
                    PbodyHydro(mycount)->U = Uinit;

#ifdef USE_CELIB //{
                    //CELibSetPrimordialMetallicity(Pbody[i]->Mass,PbodyHydro(i)->Elements);
                    //CELibSetSolarMetallicity(Pbody[i]->Mass,PbodyHydro(i)->Elements);
                    if(mode == 0){
                        CELibSetMetallicityWithSolarAbundancePattern(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements,
                            0.5*CELibGetMetalFractionForSolarChemicalComposision());
                        //CELibSetMetallicityWithSolarAbundancePattern(Pbody[i]->Mass,PbodyHydro(i)->Elements,0.02);
                        double MassLightElements = PbodyHydro(mycount)->Elements[CELibYield_H]+PbodyHydro(mycount)->Elements[CELibYield_He];
                        double Z = (Pbody[mycount]->Mass-MassLightElements)/Pbody[mycount]->Mass;
                        PbodyHydro(mycount)->Z = PbodyHydro(mycount)->ZII = PbodyHydro(mycount)->ZIa = Z;
                    } else if(mode == 1){
                        if(Pbody[mycount]->Pos[0] >0.e0){ // Right
                            CELibSetMetallicityWithSolarAbundancePattern(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements,
                                RightMetal*CELibGetMetalFractionForSolarChemicalComposision());
                        } else { //Left
                            CELibSetMetallicityWithSolarAbundancePattern(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements,
                                LeftMetal*CELibGetMetalFractionForSolarChemicalComposision());
                        }
                        double MassLightElements = PbodyHydro(mycount)->Elements[CELibYield_H]+PbodyHydro(mycount)->Elements[CELibYield_He];
                        double Z = (Pbody[mycount]->Mass-MassLightElements)/Pbody[mycount]->Mass;
                        PbodyHydro(mycount)->Z = PbodyHydro(mycount)->ZII = PbodyHydro(mycount)->ZIa = Z;
                    }
#else 
                    PbodyHydro(mycount)->Z   = 0.02;
                    PbodyHydro(mycount)->ZII = 0.02;
                    PbodyHydro(mycount)->ZIa = 0.00;
#endif // USE_CELIB //}

                    mycount ++;
                } 
            }
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    ActivateAllparticles();


#if 0
    for(int k=0;k<NProcs;k++){
        if(k == MyID){
            FILE *fp_check;
            char fname2[MaxCharactersInLine];
            sprintf(fname2,"velocity.Init.%02d.%02d",NProcs,MyID);
            FileOpen(fp_check,fname2,"w");

            for(int i=0;i<Pall.Nhydro;i++){
                fprintf(fp_check,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                    PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                    Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
            }
            fclose(fp_check);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 128;
    Pall.Npm =  8;

    Pall.TEnd = 10*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MyID == MPI_ROOT_RANK)
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
            Pall.TEnd,Pall.TEnd*Pall.UnitTime);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.0001*Pall.TEnd;

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/TC.ASCII");
    strcpy(Pall.BaseFileName,"./data/TC");
    strcpy(Pall.RestartFileName,"./data/TC.dump");

    return ;
}
#endif
#ifdef TASK_ROTATINGDISK_WITH_SINK
void InitRotatingDiskWithSink(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = 193*YEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pall.Nhydro = Number/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nsink = Pall.Nsink_t = 1;
        Pall.Ntotal = Pall.Nhydro + Pall.Nsink;
        Pall.Ntotal_t = Pall.Nhydro_t + Pall.Nsink_t;
    } else {
        Pall.Nhydro = Number/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nsink = 0;
        Pall.Nsink_t = 1;
        Pall.Ntotal = Pall.Nhydro;
        Pall.Ntotal_t = Pall.Nhydro_t + Pall.Nsink_t;
    }
    int AllocationSize = Pall.Ntotal; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        GenerateStructPsink(1);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    /* */
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pbody[AllocationSize-1]->Baryon = (void *)(Psink[0]);
        Psink[0]->Body = Pbody[AllocationSize-1];
    }


    double mgas = (0.001*MSUN_CGS/Pall.UnitMass)/Pall.Nhydro_t;
    fprintf(stderr,"%g [Msun], Nhydro %ld\n",mgas,Pall.Nhydro_t);
    //double eps = 1.e+16/Pall.UnitLength;
    double eps = SINKSINK_MERGING_DISTANCE/Pall.UnitLength;
    double cs = 0.1666*VELOCITY_KMS_CGS*GetUnitVel(); // 0.17 km/s
    // Checked.

    int count = 0;
    double Rmax = (2.56e+16/Pall.UnitLength);

    while(count < Pall.Nhydro){
    //while(count < Pall.Ntotal){
        double Pos[3];
        Pos[0] = Rmax*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        Pos[1] = Rmax*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        Pos[2] = Rmax*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);

        double r2 = NORM2(Pos);
        if(r2<SQ(Rmax)){

            Pbody[count]->Active = ON;
            Pbody[count]->Use = ON;
            Pbody[count]->Type = TypeHydro;
            Pbody[count]->GlobalID = count+Pall.Ntotal*MyID;

            Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
            Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
            Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Pos[2];
        
            double r = sqrt(SQ(Pos[0])+SQ(Pos[1]));
            double v_tan = 1.54e-12*r*Pall.UnitTime;

            Pbody[count]->Vel[0] = -v_tan*Pos[1]/r;
            Pbody[count]->Vel[1] =  v_tan*Pos[0]/r;
            Pbody[count]->Vel[2] = 0.e0;

            Pbody[count]->Eps = eps;
            Pbody[count]->Mass = mgas;
#if 1
            PbodyHydro(count)->Mass = Pbody[count]->Mass = mgas;
            PbodyHydro(count)->Z = 0.02;
#if (UseSFModelSpawn) 
            PbodyHydro(count)->SpawnMass = mgas/((double)MaxSpawnTimes);
#endif

            PbodyHydro(count)->Use = ON;
            PbodyHydro(count)->Kernel = eps;
            //PbodyHydro(count)->U = Uinit;
#endif
            count ++;
        } else {
        }
    }

    if(MPIGetMyID() == MPI_ROOT_RANK){
        int index = Pall.Ntotal-1;
        Pbody[index]->Active = ON;
        Pbody[index]->Use = ON;
        Pbody[index]->Type = TypeSink;
        Pbody[index]->GlobalID = Pall.Ntotal_t-1;

        Pbody[index]->Pos[0] = Pbody[index]->Pos[1] = Pbody[index]->Pos[2] = 0.0;
        Pbody[index]->Vel[0] = Pbody[index]->Vel[1] = Pbody[index]->Vel[2] = 0.e0;

        Pbody[index]->Eps = eps;
        Pbody[index]->Mass = 1.0*MSUN_CGS/Pall.UnitMass;

        PbodySink(index)->Use = ON;
        PbodySink(index)->ParentGlobalID = Pbody[index]->GlobalID;
        PbodySink(index)->FormationTime = 0.e0;
        PbodySink(index)->PosP[0] = PbodySink(index)->PosP[1] = PbodySink(index)->PosP[2] = 0.e0;
        PbodySink(index)->VelP[0] = PbodySink(index)->VelP[1] = PbodySink(index)->VelP[2] = 0.e0;

        PbodySink(index)->AccretionRadius = SINKHYDRO_ACCRETION_RADIUS/Pall.UnitLength;
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
        PhydroBody(i)->Vel[2] = cs*Gaussian();
    }
    
    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Active = ON;
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;
    }

#if 0
    for(int i=0;i<Pall.Ntotal;i++){
        double R = sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1]));
        // Selfgravity of gaseous disk.
        if(R > 0.e0){
            double r = sqrt(R*R+eps*eps);
            double Acc_Sink = Pall.GravConst*(1.e0/Pall.UnitMass)/CUBE(r);
            Pbody[i]->Acc[0] -= Acc_Sink*Pbody[i]->Pos[0];
            Pbody[i]->Acc[1] -= Acc_Sink*Pbody[i]->Pos[1];

            double Acc = sqrt(SQ(Pbody[i]->Acc[0])+SQ(Pbody[i]->Acc[1]));
            Pbody[i]->Vel[0] = -sqrt(Acc/R)*Pbody[i]->Pos[1];
            Pbody[i]->Vel[1] = +sqrt(Acc/R)*Pbody[i]->Pos[0];
        }
    }
#endif 

    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"./init.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]);
    }
    fclose(fp_init);

    ReConnectPointers();
    UpdateTotalNumber();

    //
    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.CS = cs;

    // Pall.Gamma = 1.0;
    // Pall.Gm1   = Pall.Gamma-1.0;
    // Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.TEnd = 1000;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
                Pall.TEnd,Pall.TEnd*Pall.UnitTime);
    }


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.0001*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/Rdisk.ASCII");
    strcpy(Pall.BaseFileName,"./data/Rdisk");
    strcpy(Pall.RestartFileName,"./data/ARdisk.dump");

    return;
}
#endif //TASK_ROTATINGDISK_WITH_SINK

void InitBreakingDam(const int MultiplicationConstant){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    // Set unit
    Pall.UnitLength = 1.e0;
    Pall.UnitTime = 1.e0;
    Pall.UnitMass = 1.e0;

    //Pall.GravConst = GetUnitGravitationalConstant();
    Pall.GravConst = 980.665; // cm/s^2

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    eprintlmpi(Pall.GravConst);
    eprintlmpi(Pall.DegreeOfFreedom);
    eprintlmpi(Pall.HeliumWeight);
    eprintlmpi(Pall.MeanMolecularWeight);
    eprintlmpi(Pall.ConvertUtoT);
    eprintlmpi(Pall.ConvertTtoU);

    eprintlmpi(Pall.ConvertDensityToCGS);
    eprintlmpi(Pall.ConvertNumberDensityToCGS);

    eprintlmpi(GetUnitVel());

    // Set Parameters
    double Lsize = 2500.0; // 25m
    double Zsize = 500.0;  // 5m
    double Mass = SQ(Lsize)*Zsize; // 1g/cm^3 * L * L * Z
    //double cs = 147700; // 1477 m/s
    double cs = 1000.0; // 1477 m/s
    // Isothermal run 

#define Nfact   (1)

    //int Nbase = Lsize*Lsize*Zsize/CUBE(100.e0);
    int Nbase = Zsize/100.e0;
    int Number = SQ(5*Nfact) * Nfact * Nbase * MultiplicationConstant;

    fprintf(stderr,"Number of particles = %d, Base particle number = %d\n",Number,Nbase);

    int count = 0;
    int mycount = 0;
    for(int i=0;i<5*Nfact*Nbase*MultiplicationConstant;i++){
    for(int j=0;j<Nfact*Nbase*MultiplicationConstant;j++){
    for(int k=0;k<5*Nfact*Nbase*MultiplicationConstant;k++){
        if(count%NProcs == MyID){
            mycount ++;
        }
        count ++;
    }
    }
    }

    int Nwater = count;

    // boundary size.
    int NboundaryX = (4*5)*Nfact*Nbase*MultiplicationConstant;
    int NboundaryY = Nfact*Nbase*MultiplicationConstant;
    int Nboundary = NboundaryX*NboundaryY;
    for(int i=0;i<NboundaryX;i++){
    for(int j=0;j<NboundaryY;j++){
        if(count%NProcs == MyID){
            mycount ++;
        }
        count ++;
    }
    }

    fprintf(stderr,"mycout = %d, Number = %d\n",mycount,Number);



    Pall.Lbox[0] = 50000.e0;
    Pall.Lboxh[0] = 0.5*Pall.Lbox[0];
    Pall.Lbox[1] = 500.e0;
    Pall.Lboxh[1] = 0.5*Pall.Lbox[1];
    Pall.Lbox[2] = 50000.e0;
    Pall.Lboxh[2] = 0.5*Pall.Lbox[2];

    Pall.BoxCenter[0] = Pall.BoxCenter[1] = Pall.BoxCenter[2] = 0.e0;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    PbodySize = mycount;
    PhydroSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double mass = Mass/(double)Nwater;
    double eps = TINY;
    double Uinit = SQ(cs)/(5.0/3.0 * (5.0/3.0-1.0));

    count = 0;
    mycount = 0;
    for(int i=0;i<5*Nfact*Nbase*MultiplicationConstant;i++){
    for(int j=0;j<Nfact*Nbase*MultiplicationConstant;j++){
    for(int k=0;k<5*Nfact*Nbase*MultiplicationConstant;k++){
        if(count%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeHydro;
            Pbody[mycount]->GlobalID = count;

            Pbody[mycount]->Pos[0] = (i+1)*Lsize/(5*Nfact*Nbase*MultiplicationConstant);
            Pbody[mycount]->Pos[1] = (j+1)*Zsize/(Nfact*Nbase*MultiplicationConstant);
            Pbody[mycount]->Pos[2] = (k+1)*Lsize/(5*Nfact*Nbase*MultiplicationConstant);

            Pbody[mycount]->Mass = mass;

            Pbody[mycount]->Vel[0] = 0.e0;
            Pbody[mycount]->Vel[1] = 0.e0;
            Pbody[mycount]->Vel[2] = 0.e0;

            Pbody[mycount]->Eps = TINY;

            PbodyHydro(mycount)->Use = ON;
            PbodyHydro(mycount)->Kernel = 2*Nfact*Nbase*MultiplicationConstant;

            mycount ++;
        }
        count ++;
    }
    }
    }

    // Add boundary particles
    for(int i=0;i<NboundaryX;i++){
    for(int j=0;j<NboundaryY;j++){
        if(count%NProcs == MyID){
            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeHydro;
            Pbody[mycount]->GlobalID = -1;

            Pbody[mycount]->Pos[0] = i*Lsize/(5*Nfact*Nbase*MultiplicationConstant);
            Pbody[mycount]->Pos[1] = j*Zsize/(Nfact*Nbase*MultiplicationConstant);
            Pbody[mycount]->Pos[2] = 0.e0;

            Pbody[mycount]->Mass = mass;

            Pbody[mycount]->Vel[0] = 0.e0;
            Pbody[mycount]->Vel[1] = 0.e0;
            Pbody[mycount]->Vel[2] = 0.e0;

            Pbody[mycount]->Eps = TINY;

            PbodyHydro(mycount)->Use = ON;
            PbodyHydro(mycount)->Kernel = 2*Nfact*Nbase*MultiplicationConstant;

            mycount ++;
        }
        count ++;
    }
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }


    // Check Total Mass
#if 0
    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"Dam.Init.%02d",MyID);
    FileOpen(fp,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
    }
    fclose(fp);
    exit(0);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 50;
    Pall.Npm = 2;
    Pall.CS = cs;

    Pall.TEnd = 10.e0;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
            Pall.TEnd,Pall.TEnd*Pall.UnitTime);

    Pall.Gamma = 1.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    strcpy(Pall.ASCIIFileName,"./data/Dam.ASCII");
    strcpy(Pall.BaseFileName,"./data/Dam");
    strcpy(Pall.RestartFileName,"./data/Dam.dump");
    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 1.e0;

    OutPutASCIIDATA();

    FileOutPutConstantInterval();
    InitializeRandomGenerator(1977);

    InitLogFiles();

    fflush(NULL);

    return;
}

void GravitationalForceFromEarth(void){

    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroActive(i)){
            PhydroAcc(i)[0] =
            PhydroAcc(i)[1] = 
            PhydroAcc(i)[2] = 0.e0;
            //if(PhydroBody(i)->GlobalID != -1)
                PhydroAcc(i)[2] -= Pall.GravConst; // gravitatinal acc. on earth.
        }
    }

    return ; 
}

void AntiForceFromBounaryRegions(const int MultiplicationConstant){

    const static double Lsize = 2500.e0; // 25m
    const static double Zsize = 500.e0;  // 5m
    //const static int Nbase = Zsize/100.e0;
    static int Nbase = 500.e0/100.e0;
    //double mean_separation = Lsize/(5.0*Nfact*Nbase*MultiplicationConstant);
    double mean_separation = Lsize/(5.0*Nfact*Nbase*MultiplicationConstant);

    fprint(mean_separation);

#define n1_for_boundary_force   (12.0)
#define n2_for_boundary_force   (0.5*n1_for_boundary_force)

    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroActive(i)){
        if(PhydroBody(i)->GlobalID != -1){
            if(PhydroPosP(i)[2] < mean_separation){ // boundary "z"
                //double r0_r = PhydroPosP(i)[2]/mean_separation;
                double r0_r = mean_separation/PhydroPosP(i)[2];
                PhydroPosP(i)[2] += -5*Pall.GravConst*Lsize*
                    (pow(r0_r,n1_for_boundary_force)-pow(r0_r,n2_for_boundary_force))/
                        fabs(PhydroPosP(i)[2]);
            }
        }
        }
    }

    for(int i=0;i<Pall.Nhydro;i++){
        if(PhydroActive(i)){
        if(PhydroBody(i)->GlobalID != -1){
            if(PhydroPosP(i)[0] < mean_separation){ // boundary "x"
                double r0_r = PhydroPosP(i)[0]/mean_separation;
                PhydroPosP(i)[0] += -5*Pall.GravConst*Lsize*
                    (pow(r0_r,n1_for_boundary_force)-pow(r0_r,n2_for_boundary_force))/
                        fabs(PhydroPosP(i)[0]);
            }
        }
        }
    }

    return ;
}


void InitUniformSphereTest(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    double r;
    int NX,NY,NZ;
    NX = NY = NZ = number;

    InitializeRandomGenerator(1977);
	//double hmean = 2.e0/cbrt((double)number);

    // count Pall.Ntotal,Pall.Nhydro

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
                    
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    fprintf(stderr,"[%02d] %d %d\n",MPIGetMyID(),mycount,count);
    memset(&Pall,0,sizeof(struct StructPall));
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    PbodySize = mycount;
    PhydroSize = mycount;
    PbodyElementsSize = mycount;
    PhydroElementsSize = mycount;

    dprintlmpi(count);
    dprintlmpi(mycount);
    dlprintlmpi(Pall.Ntotal);
    dprintlmpi(AllocationSize);

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = 1.e0/(double)count;
    double eps = 0.1*dx;
    //double eps = 0.1*(cbrt((double)1736/(double)count));
    //double eps = 0.01*(cbrt((double)1736/(double)count));
    double Uinit = 0.05;

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

				r = NORM(Pos);
				/*** make unit vector ***/
                if(r>TINY){
                    Pos[0] /= r;
                    Pos[1] /= r;
                    Pos[2] /= r;
				
                    /*** strech grid(r_new = r_old**3) ***/
                    Pos[0] *= (r)*sqrt(r);
                    Pos[1] *= (r)*sqrt(r);
                    Pos[2] *= (r)*sqrt(r);
                } else {
                    Pos[0] = 0.e0;
                    Pos[1] = 0.e0;
                    Pos[2] = 0.e0;
                }
				
				if(1.e0>NORM(Pos)){
                    if(count%NProcs == MyID){
                        Pbody[mycount]->Active = ON;
                        Pbody[mycount]->Use = ON;
                        Pbody[mycount]->Type = TypeHydro;
                        Pbody[mycount]->GlobalID = count;

                        Pbody[mycount]->Pos[0] = Pos[0];
                        Pbody[mycount]->Pos[1] = Pos[1];
                        Pbody[mycount]->Pos[2] = Pos[2];

                        Pbody[mycount]->Vel[0] = TINY;
                        Pbody[mycount]->Vel[1] = TINY;
                        Pbody[mycount]->Vel[2] = TINY;

                        Pbody[mycount]->Mass = mass;
                        Pbody[mycount]->Eps = eps;

                        Pbody[mycount]->InteractionList = 1;

                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = 0.1*dx;
                        //PbodyHydro(mycount)->Kernel = 1*dx;
                        PbodyHydro(mycount)->U = Uinit;

                        //Pall.Ntotal ++;
                        //Pall.Nhydro ++;

                        mycount ++;
                    }
					count ++;
				}
			}
		}
	}

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }


    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 5;

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    InitLogFiles();
    return;
}


void InitRandomSphereTest(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = number;

    int mycount = number/NProcs;

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MyID);
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = mycount*NProcs;
    int AllocationSize = mycount; 
    PbodySize = mycount;
    PhydroSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double dx = 1.e0/cbrt((double)number);
    double dk = dx;

    double mass = 1.e0/(double)Pall.Ntotal_t;
    double eps = 0.1*dx;
    double Uinit = 0.05;

    for(int i=0;i<mycount;i++){
        double Pos[3];
        do{
            Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
            Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        } while ( NORM2(Pos) > 1.0 );
        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        Pbody[i]->Pos[2] = Pos[2];

        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = i+mycount*MyID;


        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = mass;
        Pbody[i]->Eps = eps;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = dk;
        PbodyHydro(i)->U = Uinit;
    }

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 5;

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    InitLogFiles();
    return;
}

void InitUniformBoxTest(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int NX,NY,NZ;
    NX = NY = NZ = number;

	double dx = 2.e0/(double)(NX-1);
    double Pos[3];
    int count = 0;
    int mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

                if(count%NProcs == MyID){
                    mycount ++;
                }
                count ++;
			}
		}
	}

    memset(&Pall,0,sizeof(struct StructPall));
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = count;
    int AllocationSize = mycount; 
    PbodySize = mycount;
    PhydroSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double mass = 1.e0/(double)count;
    double eps = 0.1*dx;
    //double eps = 0.1*(cbrt((double)1736/(double)count));
    //double eps = 0.01*(cbrt((double)1736/(double)count));
    double Uinit = 0.05;

	dx = 2.e0/(double)(NX-1);
	count = 0; mycount = 0;
	for(int i=0;i<NX;i++){
		for(int j=0;j<NY;j++){
			for(int k=0;k<NZ;k++){
				Pos[0] = -1.e0+i*dx;
				Pos[1] = -1.e0+j*dx;
				Pos[2] = -1.e0+k*dx;

                if(count%NProcs == MyID){
                    Pbody[mycount]->Active = ON;
                    Pbody[mycount]->Use = ON;
                    Pbody[mycount]->Type = TypeHydro;
                    Pbody[mycount]->GlobalID = count;

                    Pbody[mycount]->Pos[0] = Pos[0];
                    Pbody[mycount]->Pos[1] = Pos[1];
                    Pbody[mycount]->Pos[2] = Pos[2];

                    Pbody[mycount]->Vel[0] = TINY;
                    Pbody[mycount]->Vel[1] = TINY;
                    Pbody[mycount]->Vel[2] = TINY;

                    Pbody[mycount]->Mass = mass;
                    Pbody[mycount]->Eps = eps;

                    PbodyHydro(mycount)->Use = ON;
                    PbodyHydro(mycount)->Kernel = 0.1*dx;
                    //PbodyHydro(mycount)->Kernel = 1*dx;
                    PbodyHydro(mycount)->U = Uinit;

                    mycount ++;
                }
                count ++;
			}
		}
	}


    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 5;

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    InitLogFiles();
    return;
}

void InitRandomBoxTestHydro(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int mycount = 0;
    for(int i=0;i<number;i++){
        if(i%NProcs == MyID){
            mycount ++;
        }
    }

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977);
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = mycount*NProcs;
    int AllocationSize = mycount;
    PbodySize = mycount;
    PhydroSize = mycount;
    PbodyElementsSize = mycount;
    PhydroElementsSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double dx = 1.e0/cbrt((double)number);
    double dk = dx;

    double mass = 1.e0/(double)Pall.Ntotal_t;
    double eps = 0.1*dx;
    double Uinit = 0.05;


    mycount = 0;
    for(int i=0;i<number;i++){
        double Pos[3];
        Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        if(i%NProcs == MyID){
            Pbody[mycount]->Pos[0] = Pos[0];
            Pbody[mycount]->Pos[1] = Pos[1];
            Pbody[mycount]->Pos[2] = Pos[2];

            Pbody[mycount]->Active = ON;
            Pbody[mycount]->Use = ON;
            Pbody[mycount]->Type = TypeHydro;
            Pbody[mycount]->GlobalID = i;


            Pbody[mycount]->Vel[0] = TINY;
            Pbody[mycount]->Vel[1] = TINY;
            Pbody[mycount]->Vel[2] = TINY;

            Pbody[mycount]->Mass = mass;
            Pbody[mycount]->Eps = eps;

            PbodyHydro(mycount)->Use = ON;
            PbodyHydro(mycount)->Kernel = dk;
            PbodyHydro(mycount)->U = Uinit;
            mycount ++;
        }
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    InitLogFiles();
    return;
}

void InitRandomBoxTest(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int mycount = number/NProcs;

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977);
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.NDM = mycount;
    Pall.Ntotal_t = Pall.NDM_t = mycount*NProcs;
    int AllocationSize = mycount; 
    PbodySize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));

    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<AllocationSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Pbody[i] = PbodyElements+i;

    double dx = 0.1/cbrt((double)number);

    double mass = 1.e0/(double)Pall.Ntotal_t;
    double eps = 0.1*dx;

    int Count = 0;
    for(int i=0;i<number;i++){
        double Pos[3];
        Pos[0] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        Pos[1] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        Pos[2] = 2.0*gsl_rng_uniform(RandomGenerator)-1.0;
        if(i%NProcs == MyID){
            Pbody[Count]->Pos[0] = Pos[0];
            Pbody[Count]->Pos[1] = Pos[1];
            Pbody[Count]->Pos[2] = Pos[2];

            Pbody[Count]->Active = ON;
            Pbody[Count]->Use = ON;
            Pbody[Count]->Type = TypeDM;
            //Pbody[Count]->GlobalID = i+mycount*MyID;
            Pbody[Count]->GlobalID = i;

            Pbody[Count]->Vel[0] = TINY;
            Pbody[Count]->Vel[1] = TINY;
            Pbody[Count]->Vel[2] = TINY;

            Pbody[Count]->Mass = mass;
            Pbody[Count]->Eps = eps;
            Pbody[Count]->InteractionList = 1;
            Count ++;
        }
    }


    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.GravConst  = 1.e0;

    InitLogFiles();
    return;
}


// For keosoku
struct StructIOKeisoku{
    long int ID;
    double Pos[3];
    double Vel[3];
    double Mass;
    double Eps;
    double KernelSize;
    double Uinit;
};

void ReadKeisokuCheckParallel(const int number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int mycount = 0;
    for(int i=0;i<number;i++){
        if(i%NProcs == MyID){
            mycount ++;
        }
    }

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977);
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = number;
    int AllocationSize = mycount;
    PbodySize = mycount;
    PhydroSize = mycount;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(AllocationSize*sizeof(StructPhydro));
    Phydro = malloc(PbodySize*sizeof(StructPhydroptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
    memset(PhydroElements,0,AllocationSize*sizeof(StructPhydro));

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize-1;i++){
        PbodyElements[i].Next = &(PbodyElements[i+1]);
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    }
    PbodyElements[AllocationSize-1].Next = NULL;
    PhydroElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++){
        Pbody[i] = PbodyElements+i;
        Phydro[i] = PhydroElements+i;
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    for(int k=0;k<NProcs;k++){
        if(k == MyID){
            FILE *fp;
            char fname[MaxCharactersInLine]; 
            sprintf(fname,"InitialCondition.Binary.%08ld.01.00",(long int)number);
            FileOpen(fp,fname,"rb");
            struct StructIOKeisoku IOTEMP;

            mycount = 0;
            for(int i=0;i<number;i++){
                long int ID;
                double Pos[3];
                double Vel[3];
                double Mass,Eps;
                double KernelSize,Uinit;
                //fscanf(fp,"%ld %le %le %le %le %le %le %le %le %le %le\n",
                        //&ID,Pos,Pos+1,Pos+2,Vel,Vel+1,Vel+2,&Eps,&Mass,&KernelSize,&Uinit);
                     
                fread(&IOTEMP,sizeof(struct StructIOKeisoku),1,fp);

                ID = IOTEMP.ID;
                Pos[0] = IOTEMP.Pos[0];
                Pos[1] = IOTEMP.Pos[1];
                Pos[2] = IOTEMP.Pos[2];

                Vel[0] = IOTEMP.Vel[0];
                Vel[1] = IOTEMP.Vel[1];
                Vel[2] = IOTEMP.Vel[2];

                Eps = IOTEMP.Eps;
                Mass = IOTEMP.Mass;
                KernelSize = IOTEMP.KernelSize;
                Uinit = IOTEMP.Uinit;

                if(i%NProcs == MyID){
                    Pbody[mycount]->Pos[0] = Pos[0];
                    Pbody[mycount]->Pos[1] = Pos[1];
                    Pbody[mycount]->Pos[2] = Pos[2];

                    Pbody[mycount]->Vel[0] = Vel[0];
                    Pbody[mycount]->Vel[1] = Vel[1];
                    Pbody[mycount]->Vel[2] = Vel[2];

                    Pbody[mycount]->Active = ON;
                    Pbody[mycount]->Use = ON;
                    Pbody[mycount]->Type = TypeHydro;
                    Pbody[mycount]->GlobalID = ID;


                    Pbody[mycount]->Mass = Mass;
                    Pbody[mycount]->Eps = Eps;

                    PbodyHydro(mycount)->Use = ON;
                    PbodyHydro(mycount)->Kernel = KernelSize;
                    PbodyHydro(mycount)->U = Uinit;
                    mycount ++;
                }
            }
            fclose(fp);
        }
    }

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"Init.%02d.%02d",NProcs,MyID);
    FileOpen(fp ,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                PbodyHydro(i)->Kernel);
    }
    fclose(fp);
#endif
    
    Pall.RunStatus = NewSimulation;

    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = 3.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst  = 1.e0;

    return;
}


void ReadKeisokuBenchParallel(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitLength = MPC_CGS;
    Pall.UnitTime = 10.0*GIGAYEAR_CGS;
    Pall.UnitMass = 1.e+11*MSUN_CGS;
    Pall.TCMB = CMB_TEMPERATURE;


    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 10.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    for(int k=0;k<NProcs;k++){
        if(k == MyID){
            FILE *fp;
            char fname[MaxCharactersInLine];
            sprintf(fname,"InitialCondition.Binary.Cosmological.02011200.01.00");
            FileOpen(fp,fname,"rb");
            fread(&Pall.Ntotal_t,sizeof( unsigned long int),1,fp);
            fread(&Pall.Nhydro_t,sizeof( unsigned long int),1,fp);
            fclose(fp);

            Pall.NDM_t = Pall.Ntotal_t-Pall.Nhydro_t;

            Pall.NDM = 0;
            Pall.Nhydro = 0;
            Pall.Ntotal = 0;
            for(int i=0;i<Pall.Ntotal_t;i++){
                if(i%NProcs == MyID){
                    if(i < Pall.Nhydro_t)
                        Pall.Nhydro ++;
                    else
                        Pall.NDM ++;
                    Pall.Ntotal ++;
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }



    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    PbodySize = Pall.Ntotal;
    PhydroSize = Pall.Nhydro;

    dlprintlmpi(Pall.NDM);
    dlprintlmpi(Pall.Ntotal);
    dlprintlmpi(Pall.Nhydro);

    PbodyElements = malloc(PbodySize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    PhydroElements = malloc(PhydroSize*sizeof(StructPhydro));
    Phydro = malloc(PhydroSize*sizeof(StructPhydroptr));

    memset(PbodyElements,0,PbodySize*sizeof(StructPbody));
    memset(PhydroElements,0,PhydroSize*sizeof(StructPhydro));


    for(int i=0;i<PbodySize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[PbodySize-1].Next = NULL;

    for(int i=0;i<PbodySize;i++)
        Pbody[i] = PbodyElements+i;


    for(int i=0;i<PhydroSize-1;i++)
        PhydroElements[i].Next = &(PhydroElements[i+1]);
    PhydroElements[PhydroSize-1].Next = NULL;

    for(int i=0;i<PhydroSize;i++)
        Phydro[i] = PhydroElements+i;

    for(int i=0;i<PhydroSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    for(int k=0;k<NProcs;k++){
        if(k == MyID){
            FILE *fp;
            char fname[MaxCharactersInLine];
            sprintf(fname,"InitialCondition.Binary.Cosmological.02011200.01.00");
            FileOpen(fp,fname,"rb");
            unsigned long int idummy;
            fread(&idummy,sizeof( unsigned long int),1,fp);
            fread(&idummy,sizeof( unsigned long int),1,fp);

            struct StructIOKeisoku IOTEMP;
            int mycount = 0;
            for(int i=0;i<Pall.Ntotal_t;i++){
                fread(&IOTEMP,sizeof(struct StructIOKeisoku),1,fp);
                if(i%NProcs == MyID){
                    Pbody[mycount]->Pos[0] = IOTEMP.Pos[0];
                    Pbody[mycount]->Pos[1] = IOTEMP.Pos[1];
                    Pbody[mycount]->Pos[2] = IOTEMP.Pos[2];

                    Pbody[mycount]->Vel[0] = IOTEMP.Vel[0];
                    Pbody[mycount]->Vel[1] = IOTEMP.Vel[1];
                    Pbody[mycount]->Vel[2] = IOTEMP.Vel[2];

                    Pbody[mycount]->Mass = IOTEMP.Mass;
                    Pbody[mycount]->Eps = IOTEMP.Eps;

                    Pbody[mycount]->GlobalID = IOTEMP.ID;

                    Pbody[mycount]->Active = ON;
                    Pbody[mycount]->Use = ON;

                    if(i<Pall.Nhydro_t){
                        PbodyHydro(mycount)->Use = ON;
                        PbodyHydro(mycount)->Kernel = IOTEMP.KernelSize;
                        PbodyHydro(mycount)->U = IOTEMP.Uinit;
                        Pbody[mycount]->Type = TypeHydro;
                    }else{
                        Pbody[mycount]->Type = TypeDM;
                    }
                    mycount ++;
                }
            }
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"Init.%02d.%02d",NProcs,MyID);
    FileOpen(fp ,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                PbodyHydro(i)->Kernel);
    }
    fclose(fp);
#endif
    Pall.RunStatus = NewSimulation;

    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.hubble = 0.01*50;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    Pall.OmegaM = 1.0;
    Pall.OmegaL = 0.0; 

    Pall.TEnd = CalcZtoT(0.e0);
    Pall.TCurrent = 0.481333;
    Pall.Redshift = CalcTtoZ(Pall.TCurrent);
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    UpdateAdaptiveSofteningFactor();

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
        Pall.AdaptiveSofteningFactor);

    //InitLogFiles();

    return;
}

#ifdef TASK_SANTABARBARA

static int CountTagZero(FILE *fp, const int Number){

    int counter;
    for(int i=0;i<Number;i++){
        int Tag;
        fscanf(fp,"%*e %*e %*e %*e %*e %*e %*e %d",&Tag);
        if(Tag == 0){
            counter ++;
        }
    }

    return counter;
}

/*
 * This function reads the initial condition file of the Santa Barabara cluster
 * comparision test. The parameter, mode, is the switch for controlling the run mode.
 * mode == 0  is an N-body only run, whereas mode == 1 is an N-body/SPH run.
 */
void ReadSantaBarbaraInitialCondition(const int mode){ 

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = MPC_CGS;
    Pall.UnitTime = (MPC_CGS/1.e+5); // Velocity unit is set to 1km/s
    Pall.UnitMass = MSUN_CGS;
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    // Set cosmological parameters!

    FILE *fp;
    int Number;
    float dummy[5];
    FileOpen(fp,"SantaBarbara.dat","r");

    fscanf(fp,"%e %e %e %e %e",dummy,dummy+1,dummy+2,dummy+3,dummy+4);
    fscanf(fp,"%d",&Number);
    if(MyID == MPI_ROOT_RANK){
        fprintf(stderr,"H0 = %g\n",dummy[0]);
        fprintf(stderr,"Omega_m = %g\n",dummy[1]);
        fprintf(stderr,"Omega_b = %g\n",dummy[2]);
        fprintf(stderr,"Omega_l = %g\n",dummy[3]);
        fprintf(stderr,"Lbox = %g\n",dummy[4]);
        fprintf(stderr,"Number of particles = %d\n",Number);
        //dprintlmpi(Number);
    }

    bool MultiMassInitCondition = false;
    int NumberTagZero = 0;
    if(Number < 0) {
        MultiMassInitCondition = true;
        Number =  -Number;
        // fpos_t *pos;
        // int fp_log = fgetpos(fp,pos);
        // NumberTagZero = CountTagZero(fp,Number);
        // fp_log = fsetpos(fp,pos);
    }

    // Allocate part
    int MyNumber = Number/NProcs;
    if(MyID == MPI_ROOT_RANK){
        MyNumber += Number%NProcs;
    }
    /*
    MyNumber = 0;
    for(int i=0;i<Number;i++){
        if(i%NProcs == MyID)
            MyNumber ++;
    }
    */
    dprintlmpi(MyNumber);

    if(mode == 0){
        Pall.NDM = MyNumber;
        Pall.Ntotal = Pall.NDM;
        Pall.Ntotal_t = Pall.NDM_t = Number;

        GenerateStructPbody(Pall.Ntotal);
    } else if (mode == 1){
        Pall.Nhydro = MyNumber;
        Pall.NDM = MyNumber;
        Pall.Ntotal = Pall.Nhydro+Pall.NDM;
        Pall.Nhydro_t = Number;
        Pall.NDM_t = Number;
        Pall.Ntotal_t = Pall.Nhydro_t+Pall.Ntotal_t;
        GenerateStructPbody(Pall.Ntotal);
        GenerateStructPhydro(Pall.Nhydro);
        for(int i=0;i<Pall.Nhydro;i++){
            Pbody[i]->Baryon = (void *)(Phydro[i]);
            Phydro[i]->Body = Pbody[i];
        }
    } else {
        fprintf(stderr,"Mode %d has a wrong value.\n",mode);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    InitializeAllActiveParticleNumbers();

    // double EPS_CDM = 20*KPC_CGS/Pall.UnitLength;
    // double EPS_Hydro = 10*KPC_CGS/Pall.UnitLength;
    double EPS_CDM = 10*KPC_CGS/Pall.UnitLength;
    double EPS_Hydro = 5*KPC_CGS/Pall.UnitLength;
    // double EPS_CDM = 50*KPC_CGS/Pall.UnitLength;
    // double EPS_Hydro = 25*KPC_CGS/Pall.UnitLength;
    
    if(mode == 0){
        int InitNumber = MyID*(Number/NProcs);
        if(MyID != 0) InitNumber += Number%NProcs;
        
        // fprintf(stderr,"%d -- %d %d | %d / %d\n",MPIGetMyID(),InitNumber,InitNumber+MyNumber,MyNumber,Number);
        int counter = 0;
        for(int i=0;i<Number;i++){
            float DataStream[7];
            int Tag;
            int EpsFactor = 1;
            if(MultiMassInitCondition == true){
                fscanf(fp,"%e %e %e %e %e %e %e %d",DataStream,
                        DataStream+1,DataStream+2,DataStream+3,
                        DataStream+4,DataStream+5,DataStream+6,&Tag);
                EpsFactor = 1<<(Tag+1);
            } else {
                fscanf(fp,"%e %e %e %e %e %e %e",DataStream,
                        DataStream+1,DataStream+2,DataStream+3,
                        DataStream+4,DataStream+5,DataStream+6);
            }
            if((i>=InitNumber)&&(i<InitNumber+MyNumber)){

                Pbody[counter]->Pos[0] = DataStream[0];
                Pbody[counter]->Pos[1] = DataStream[1];
                Pbody[counter]->Pos[2] = DataStream[2];
                Pbody[counter]->Vel[0] = DataStream[3];
                Pbody[counter]->Vel[1] = DataStream[4];
                Pbody[counter]->Vel[2] = DataStream[5];
                Pbody[counter]->Mass = DataStream[6];

                /*
                fprintf(stderr,"%g %g %g %g %g %g\n",
                        Pbody[counter]->Pos[0],Pbody[counter]->Pos[1],Pbody[counter]->Pos[2],
                        Pbody[counter]->Vel[0],Pbody[counter]->Vel[1],Pbody[counter]->Vel[2]);
                */

                Pbody[counter]->Active = ON;
                Pbody[counter]->Use = ON;
                Pbody[counter]->Type = TypeDM;
                Pbody[counter]->GlobalID = i;

                Pbody[counter]->Eps = EpsFactor*EPS_CDM;

                counter ++;
            }
        }
    } else if(mode == 1){ 

        int InitNumber = MyID*(Number/NProcs);
        if(MyID != 0) InitNumber += Number%NProcs;

        double Uinit = Pall.TCMB*SQ(SQ(1.0+Pall.Redshift));

        int counter = 0;
        for(int i=0;i<Number;i++){
            float DataStream[7];
            int Tag;
            int EpsFactor = 1;
            if(MultiMassInitCondition == true){
                fscanf(fp,"%e %e %e %e %e %e %e %d",DataStream,
                        DataStream+1,DataStream+2,DataStream+3,
                        DataStream+4,DataStream+5,DataStream+6,&Tag);
                EpsFactor = 1<<(Tag+1);
            } else {
                fscanf(fp,"%e %e %e %e %e %e %e",DataStream,
                        DataStream+1,DataStream+2,DataStream+3,
                        DataStream+4,DataStream+5,DataStream+6);
            }
            if((i>=InitNumber)&&(i<InitNumber+MyNumber)){

                // Hydro
                Pbody[counter]->Pos[0] = DataStream[0];
                Pbody[counter]->Pos[1] = DataStream[1];
                Pbody[counter]->Pos[2] = DataStream[2];
                Pbody[counter]->Vel[0] = DataStream[3];
                Pbody[counter]->Vel[1] = DataStream[4];
                Pbody[counter]->Vel[2] = DataStream[5];
                Pbody[counter]->Mass = 0.1*DataStream[6];

                Pbody[counter]->Active = ON;
                Pbody[counter]->Use = ON;
                Pbody[counter]->Type = TypeHydro;
                Pbody[counter]->GlobalID = i;

                Pbody[counter]->Eps = EpsFactor*EPS_Hydro;

                PbodyHydro(counter)->Mass = Pbody[counter]->Mass;
                PbodyHydro(counter)->Z = 0.0;
                PbodyHydro(counter)->ZII = 
                PbodyHydro(counter)->ZIa = 0.0;
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
                PbodyHydro(counter)->G0 = 0.0;
                PbodyHydro(counter)->fH2 = 0.1;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

                PbodyHydro(counter)->DQheat = 0.e0;
                PbodyHydro(counter)->Use = ON;
                PbodyHydro(counter)->Kernel = EpsFactor*EPS_Hydro;
                PbodyHydro(counter)->U = Uinit;
#ifdef USE_PARTICLE_TAG
                PbodyHydro(counter)->Tag = Tag;
#endif // USE_PARTICLE_TAG


                // DM
                Pbody[counter+MyNumber]->Pos[0] = DataStream[0];
                Pbody[counter+MyNumber]->Pos[1] = DataStream[1];
                Pbody[counter+MyNumber]->Pos[2] = DataStream[2];
                Pbody[counter+MyNumber]->Vel[0] = DataStream[3];
                Pbody[counter+MyNumber]->Vel[1] = DataStream[4];
                Pbody[counter+MyNumber]->Vel[2] = DataStream[5];
                Pbody[counter+MyNumber]->Mass = 0.9*DataStream[6];

                Pbody[counter+MyNumber]->Active = ON;
                Pbody[counter+MyNumber]->Use = ON;
                Pbody[counter+MyNumber]->Type = TypeDM;
                Pbody[counter+MyNumber]->GlobalID = i+Number;

                Pbody[counter+MyNumber]->Eps = EpsFactor*EPS_CDM;

                counter ++;
            }
        }
    }

    
    ActivateAllparticles();
    // UpdateTotalNumber();
    // ReConnectPointers();
    // UpdateTotalActiveNumber();

    //write cloud.
    char fcloud[MaxCharactersInLine];
    Snprintf(fcloud,"_hydro.%03d",MyID);
    FileOpen(fp,fcloud,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%d %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2]);
    }
    fclose(fp);
    fflush(NULL);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat _hydro.* > hydro.dat");
        fflush(NULL);
        system("rm -rf ./_hydro.*");
        fflush(NULL);
    }
    fflush(NULL);

    Snprintf(fcloud,"_DM.%03d",MyID);
    FileOpen(fp,fcloud,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            fprintf(fp,"%d %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
        }
    }
    fclose(fp);
    fflush(NULL);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat _DM.* > DM.dat");
        fflush(NULL);
        system("rm -rf ./_DM.*");
        fflush(NULL);
    }
    fflush(NULL);

    MPI_Barrier(MPI_COMM_WORLD);

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    //Pall.Ns = 32;
    // Pall.Ns = 200;
    // Pall.Npm = 10;
    Pall.Ns = 128;
    Pall.Npm = 8;
    //Pall.Ns = 64;
    //Pall.Npm = 5;
    // Pall.Ns = 32;
    // Pall.Npm = 2;


    Pall.hubble = 0.5;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    if(mode == 0){
        Pall.OmegaM = 1.0;
        Pall.OmegaB = 0.0;
    } else if(mode == 1){
        //Pall.OmegaM = 0.9;
        Pall.OmegaM = 1.0;
        Pall.OmegaB = 0.1;
    }
    Pall.OmegaL = 0.0;
    Pall.InitialRedshift = Pall.Redshift = 63.0;

    Pall.TEnd = CalcZtoT(0.e0);
    Pall.TCurrent = CalcZtoT(Pall.InitialRedshift);
    Pall.Redshift = Pall.InitialRedshift;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
            Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);
    UpdateAdaptiveSofteningFactor();

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.01*Pall.TEnd;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(!CheckDir("./data"))
            MakeDir("./data");
    }

    strcpy(Pall.ASCIIFileName,"./data/SantaBarbara.ASCII");
    strcpy(Pall.BaseFileName,"./data/SantaBarbara");
    strcpy(Pall.RestartFileName,"./data/SantaBarbara.dump");

    return;
}

#endif // TASK_SANTABARBARA


#ifdef TASK_GALAXY_FORMATION //{

#define USE_EPS_FOLLOW_MASS
#define USE_BARYON_SOFTENING

static double ComputeEps(const double Eps, const double Mass, const int Type){

#ifdef USE_EPS_FOLLOW_MASS //{

#define Mvir (1.e+12)
#ifdef USE_BARYON_SOFTENING //{
    if(Type == TypeHydro){
        // 100*mass/eps^3 = rhoth;
        // eps = (100*mass/rhoth)^{1/3}
        // need to consider 4pi/3 ?
        double rhoth = SFCONDITION_DENSITY_CRITERION/Pall.ConvertNumberDensityToCGS;
        return cbrt(100*Mass/rhoth);
    } else {
        double Mass_in_Solar = Mass*Pall.UnitMass/MSUN_CGS;
        double Eps_new = 30*sqrt(Mass_in_Solar/1.e+3)*pow((Mvir/1.e+12),0.2)*PC_CGS/Pall.UnitLength;
        return Eps_new/3.0;
    }
#else  // USE_BARYON_SOFTENING //}//{
    double Mass_in_Solar = Mass*Pall.UnitMass/MSUN_CGS;
    double Eps_new = 30*sqrt(Mass_in_Solar/1.e+3)*pow((Mvir/1.e+12),0.2)*PC_CGS/Pall.UnitLength;
    return Eps_new/3.0;
#endif // USE_BARYON_SOFTENING //}


#else // USE_EPS_FOLLOW_MASS //}//{
    return Eps;
#endif  // USE_EPS_FOLLOW_MASS //}
#undef Mvir

}


static void InsertData(const int TargetIndex, float Pos[], float Vel[], float Mass, float Rho, float U,
        const int Type, const int GlobalID, const double Eps){

    Pbody[TargetIndex]->GlobalID = GlobalID;
    Pbody[TargetIndex]->Pos[0] = Pos[0];
    Pbody[TargetIndex]->Pos[1] = Pos[1];
    Pbody[TargetIndex]->Pos[2] = Pos[2];

    Pbody[TargetIndex]->Vel[0] = Vel[0];
    Pbody[TargetIndex]->Vel[1] = Vel[1];
    Pbody[TargetIndex]->Vel[2] = Vel[2];
    
    Pbody[TargetIndex]->Mass = Mass;

    Pbody[TargetIndex]->Active = ON;
    Pbody[TargetIndex]->Use = ON;
    Pbody[TargetIndex]->Type = Type;
    //Pbody[TargetIndex]->Eps = Eps;
    Pbody[TargetIndex]->Eps = ComputeEps(Eps,Mass,Type);

    if(Type == TypeHydro){
        PbodyHydro(TargetIndex)->Use = ON;
        PbodyHydro(TargetIndex)->Rho = Rho;
        PbodyHydro(TargetIndex)->U = U;
        PbodyHydro(TargetIndex)->Kernel = 2.0*Eps;
#if (UseSFModelSpawn) 
        PbodyHydro(TargetIndex)->SpawnMass = Pbody[TargetIndex]->Mass/(double)MaxSpawnTimes;
#endif
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
        PbodyHydro(TargetIndex)->G0 = 1;
        PbodyHydro(TargetIndex)->fH2 = 0.1;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

#ifdef USE_CELIB //{
        CELibSetPrimordialMetallicity(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements);
        //CELibSetSolarMetallicity(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements);
        //CELibSetMetallicityWithSolarAbundancePattern(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements,
                // 0.1*CELibGetMetalFractionForSolarChemicalComposision());
        double MassLightElements = PbodyHydro(TargetIndex)->Elements[CELibYield_H]
                                  +PbodyHydro(TargetIndex)->Elements[CELibYield_He];
        double Z = (Pbody[TargetIndex]->Mass-MassLightElements)/Pbody[TargetIndex]->Mass;
        PbodyHydro(TargetIndex)->Z = PbodyHydro(TargetIndex)->ZII = PbodyHydro(TargetIndex)->ZIa = Z;
#else 
        PbodyHydro(TargetIndex)->Z   = 0.0;
        PbodyHydro(TargetIndex)->ZII = 0.0;
        PbodyHydro(TargetIndex)->ZIa = 0.0;
#endif // USE_CELIB //}
    }

    return ;
}

static void ReadMusicLine(char fname[], const int Nlines, const int mycount_comp[], const double Eps[], const double Mdm){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int Number = mycount_comp[0] + mycount_comp[1] + mycount_comp[2];
    // dprintlmpi(Number);
    // dprintlmpi(mycount_comp[0]);
    // dprintlmpi(mycount_comp[1]);
    // dprintlmpi(mycount_comp[2]);

    FILE *fp;
    FileOpen(fp,fname,"r");
    fscanf(fp,"%*d");
    fscanf(fp,"%*d");
    fscanf(fp,"%*d");
    fscanf(fp,"%*e");
    fscanf(fp,"%*e");
    fscanf(fp,"%*e");
    fscanf(fp,"%*e");
    fscanf(fp,"%*e");
    fscanf(fp,"%*e");
    fscanf(fp,"%*e");
    fscanf(fp,"%*e");

    double InitT = Pall.TCMB*SQ(SQ(1.0+Pall.InitialRedshift));
    double InitU = Pall.ConvertTtoU*InitT;

    int counter_total = 0;
    int counter_hydro = 0;
    int counter_dm = 0;
    int counter_boundary = 0;
    for(int i=0;i<Nlines;i++){
        float Pos[3],Vel[3],Mass;
        int Type;
        float Rho = 0.e0;
        fscanf(fp,"%e %e %e %e %e %e %e %d",
                Pos,Pos+1,Pos+2,Vel,Vel+1,Vel+2,
                &Mass,&Type);
#if 1
        Vel[0] = Vel[0]*1.e5/(Pall.UnitLength/Pall.UnitTime); // km/s -> cm/s -> L/T
        Vel[1] = Vel[1]*1.e5/(Pall.UnitLength/Pall.UnitTime);
        Vel[2] = Vel[2]*1.e5/(Pall.UnitLength/Pall.UnitTime);
#endif
        if(i%NProcs == MyID){
            if(Type == 0){ // hydro
                InsertData(counter_hydro,Pos,Vel,Mass,Rho,InitU,TypeHydro,i,Eps[Type]);
                counter_hydro ++;
            } else if(Type == 1){ //DM
                InsertData(counter_dm+mycount_comp[0],Pos,Vel,Mass,0.e0,0.e0,TypeDM,i,Eps[Type]);
                counter_dm ++;
            } else {
                // softening.
                double eps = 4*Eps[1]*cbrt(Mass/Mdm);
                InsertData(counter_boundary+mycount_comp[0]+mycount_comp[1],Pos,Vel,Mass,0.e0,0.e0,TypeDM,i,eps);
                counter_boundary ++;
            }
            counter_total ++;
        }
    }
    fclose(fp);

    // dprintlmpi(counter_hydro);
    // dprintlmpi(counter_dm);
    // dprintlmpi(counter_boundary);
    // dprintlmpi(Number);
    assert(Number == counter_total);
    // exit(1);

    return ;
}

static void ReadMusicLineBinary(char fname[], const int Nlines, const int mycount_comp[], const double Eps[], const double Mdm){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int Number = mycount_comp[0] + mycount_comp[1] + mycount_comp[2];

    FILE *fp;
    FileOpen(fp,fname,"rb");
    int IntDummy[3];
    double DoubleDummy[8];
    fread(IntDummy,sizeof(int),3,fp);
    fread(DoubleDummy,sizeof(double),8,fp);

    double InitT = Pall.TCMB*SQ(SQ(1.0+Pall.InitialRedshift));
    double InitU = Pall.ConvertTtoU*InitT;

    int counter_total = 0;
    int counter_hydro = 0;
    int counter_dm = 0;
    int counter_boundary = 0;
    for(int i=0;i<Nlines;i++){
        // float Pos[3],Vel[3],Mass;
        // int Type;
        float Rho = 0.e0;

        struct StructTmpWrite{
            float Pos[3];
            float Vel[3];
            float Mass;
            int Type;
        } TmpRead;

        fread(&TmpRead,sizeof(struct StructTmpWrite),1,fp);

        /*
        fprintf(stderr,"%g %g %g %g %g %g %g %d\n",
                TmpRead.Pos[0],
                TmpRead.Pos[1],
                TmpRead.Pos[2],
                TmpRead.Vel[0],
                TmpRead.Vel[1],
                TmpRead.Vel[2],
                TmpRead.Mass,TmpRead.Type);
                */
        // fscanf(fp,"%e %e %e %e %e %e %e %d",
                // Pos,Pos+1,Pos+2,Vel,Vel+1,Vel+2,
                // &Mass,&Type);
#if 1
        TmpRead.Vel[0] = TmpRead.Vel[0]*1.e5/(Pall.UnitLength/Pall.UnitTime); // km/s -> cm/s -> L/T
        TmpRead.Vel[1] = TmpRead.Vel[1]*1.e5/(Pall.UnitLength/Pall.UnitTime);
        TmpRead.Vel[2] = TmpRead.Vel[2]*1.e5/(Pall.UnitLength/Pall.UnitTime);
#endif
        if(i%NProcs == MyID){
            if(TmpRead.Type == 0){ // hydro
                InsertData(counter_hydro,TmpRead.Pos,TmpRead.Vel,TmpRead.Mass,Rho,InitU,TypeHydro,i,Eps[TmpRead.Type]);
                counter_hydro ++;
            } else if(TmpRead.Type == 1){ //DM
                InsertData(counter_dm+mycount_comp[0],TmpRead.Pos,TmpRead.Vel,TmpRead.Mass,0.e0,0.e0,TypeDM,i,Eps[TmpRead.Type]);
                counter_dm ++;
            } else {
                // softening.
                double eps = 4*Eps[1]*cbrt(TmpRead.Mass/Mdm);
                InsertData(counter_boundary+mycount_comp[0]+mycount_comp[1],TmpRead.Pos,TmpRead.Vel,TmpRead.Mass,0.e0,0.e0,TypeDM,i,eps);
                counter_boundary ++;
            }
            counter_total ++;
        }
    }
    fclose(fp);

    // dprintlmpi(counter_hydro);
    // dprintlmpi(counter_dm);
    // dprintlmpi(counter_boundary);
    // dprintlmpi(Number);
    assert(Number == counter_total);
    // exit(1);

    return ;
}

void NumberCounter(int Numbers[], char fname[]){

    for(int i=0;i<8;i++)
        Numbers[i] = 0;

    FILE *fp;
    FileOpen(fp,fname,"r");
    int Lines[4];
    fscanf(fp,"%d %d %d %d",Lines,Lines+1,Lines+2,Lines+3);
    int Nlines = Lines[0]+Lines[1]+Lines[2]+Lines[3];
    fscanf(fp,"%*e %*e %*e %*e");
    int counter = 0;
    int OldType = 0;
    for(int i=0;i<Nlines;i++){
        int Type;
        fscanf(fp,"%*e %*e %*e %*e %*e %*e %*e %*e %*e %d %*e %*e",&Type);
        if(i>0){
            if(OldType != Type)
                counter ++;
            Numbers[counter] ++;
        } else {
            Numbers[counter] ++;
        }
        OldType = Type;
    }
    fclose(fp);
    for(int i=0;i<8;i++)
        dprintlmpi(Numbers[i]);

    fflush(NULL);
    exit(1);
    return ;
}

static void CheckEpsValues(void){

    double eps_hydro = 0;
    double mass_hydro = 0;

    if(Pall.Nhydro > 0){
        eps_hydro = PhydroBody(0)->Eps;
        mass_hydro = PhydroBody(0)->Mass;
    }
    double eps_hydro_global;
    double mass_hydro_global;
    MPI_Reduce(&eps_hydro,&eps_hydro_global,1,MPI_DOUBLE,MPI_MAX,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Reduce(&mass_hydro,&mass_hydro_global,1,MPI_DOUBLE,MPI_MAX,MPI_ROOT_RANK,MPI_COMM_WORLD);


    double eps_DM[2] = {0.e0,0.e0};
    double mass_DM[2] = {0.e0,0.e0};
    int counter = 0;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            if(counter == 0){
                eps_DM[0] = eps_DM[1] = Pbody[i]->Eps;
                mass_DM[0] = mass_DM[1] = Pbody[i]->Mass;
            } else {
                eps_DM[0] = fmin(Pbody[i]->Eps,eps_DM[0]);
                eps_DM[1] = fmax(Pbody[i]->Eps,eps_DM[1]);
                mass_DM[0] = fmin(Pbody[i]->Eps,mass_DM[0]);
                mass_DM[1] = fmax(Pbody[i]->Eps,mass_DM[1]);
            }
            counter++;
        }
    }

    eps_DM[0] *= -1.e0;
    mass_DM[0] *= -1.e0;

    double eps_DM_global[2];
    double mass_DM_global[2];
    MPI_Reduce(&eps_DM[0],&eps_DM_global[0],1,MPI_DOUBLE,MPI_MIN,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Reduce(&eps_DM[1],&eps_DM_global[1],1,MPI_DOUBLE,MPI_MAX,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Reduce(&mass_DM[0],&mass_DM_global[0],1,MPI_DOUBLE,MPI_MIN,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Reduce(&mass_DM[1],&mass_DM_global[1],1,MPI_DOUBLE,MPI_MAX,MPI_ROOT_RANK,MPI_COMM_WORLD);


    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Softening report\n");
        fprintf(stderr,"Gas %g pc (%g Msun)\n",
                eps_hydro_global*Pall.UnitLength/PC_CGS,
                mass_hydro_global*Pall.UnitMass/MSUN_CGS);
        fprintf(stderr,"DM from %g pc (%g Msun) to %g pc (%g Msun)\n",
                -eps_DM_global[0]*Pall.UnitLength/PC_CGS,
                -mass_DM_global[0]*Pall.UnitMass/MSUN_CGS,
                eps_DM_global[1]*Pall.UnitLength/PC_CGS,
                mass_DM_global[1]*Pall.UnitMass/MSUN_CGS);
    }

    return ;
}

void ReadMusic(char fname[], const bool BinaryFlag){

// #ifdef USE_CELIB
    // InitializeStellarFeedback();
// #endif

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    
    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = MPC_CGS;
    //Pall.UnitTime = (MPC_CGS/1.e+5); // Velocity unit is set to 1km/s
    Pall.UnitTime = MPC_CGS/(VELOCITY_KMS_CGS);
    Pall.UnitMass = 1.e+10*MSUN_CGS;
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    //SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
    SetVariableViscosityParameters(VARIABLE_VISCOSITY_MIN,VARIABLE_VISCOSITY_MAX,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    int Numbers[3] = {0,0,0};
    double Params[8];
    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(BinaryFlag){
            FILE *fp;
            FileOpen(fp,fname,"rb");
            fread(Numbers,sizeof(int),3,fp);
            fread(Params,sizeof(double),8,fp);
            fclose(fp);
        } else {
            FILE *fp;
            FileOpen(fp,fname,"r");
            fscanf(fp,"%d",Numbers);
            fscanf(fp,"%d",Numbers+1);
            fscanf(fp,"%d",Numbers+2);
            // fscanf(fp,"%le",Time);
            // fscanf(fp,"%le",Redshift);
            // fscanf(fp,"%le",BoxSize);
            // fscanf(fp,"%le",Omega0);
            // fscanf(fp,"%le",OmegaL);
            // fscanf(fp,"%le",Hubble);
            fscanf(fp,"%le",Params);
            fscanf(fp,"%le",Params+1);
            fscanf(fp,"%le",Params+2);
            fscanf(fp,"%le",Params+3);
            fscanf(fp,"%le",Params+4);
            fscanf(fp,"%le",Params+5);
            fscanf(fp,"%le",Params+6);
            fscanf(fp,"%le",Params+7);
            fclose(fp);
        }
    }
    MPI_Bcast(&Numbers,3,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&Params,8,MPI_DOUBLE,MPI_ROOT_RANK,MPI_COMM_WORLD);

    // Set Shuffule/Lock information
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Particle numbers, Ngas,Ndm, Nboundary = %d %d %d\n",
                Numbers[0],Numbers[1],Numbers[2]);
        fprintf(stderr,"Parameters Time = %g, Redshift = %g, BoxSize = %g, Omega0 = %g, OmegaL = %g, Hubble = %g\n",
                Params[0],Params[1],Params[2],
                Params[3],Params[4],Params[5]);
        fprintf(stderr,"Parameters m_baryon = %g, m_dm = %g\n",
                Params[6],Params[7]);
    }
    Pall.Ntotal_t = Numbers[0]+Numbers[1]+Numbers[2];


    int mycount_comp[3] = {0,0,0};
    int mycount = 0;

    for(int k=0;k<Pall.Ntotal_t;k++){
        if(k%NProcs == MyID){
            if(k<Numbers[0]){
                mycount_comp[0] ++;
            } else if(k<Numbers[0]+Numbers[1]){
                mycount_comp[1] ++;
            } else if(k<Numbers[0]+Numbers[1]+Numbers[2]){
                mycount_comp[2] ++;
            } 
            mycount ++;
        }
    }


    Pall.Ntotal = mycount;
    Pall.Nhydro = mycount_comp[0];
    Pall.NDM = mycount_comp[1]+mycount_comp[2];
    Pall.Nstars = 0;

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
        Pbody[i]->Type = TypeHydro;
    }


    double Eps[] = {20.0*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {50.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength,200.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {10.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength};

    // double Eps[3];
    // Eps[TypeDM] = 100.0*PC_CGS/Pall.UnitLength; // for DM
    // Eps[TypeHydro] = 50.0*PC_CGS/Pall.UnitLength; // for hydro
    // Eps[2] = 200.0*PC_CGS/Pall.UnitLength; // for boundary 
    
#if 0
    for(int i=0;i<NProcs;i++){
        if(MyID == i){
            fprintf(stderr,"Start read file. Node ID = %d\n",i);
            if(BinaryFlag){
                ReadMusicLineBinary(fname,Numbers[0]+Numbers[1]+Numbers[2],mycount_comp,Eps,Params[7]);
            } else {
                ReadMusicLine(fname,Numbers[0]+Numbers[1]+Numbers[2],mycount_comp,Eps,Params[7]);
            }
            fprintf(stderr,"Finish read file. Node ID = %d\n",i);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#else
    for(int i=0;i<NProcs;i+=INITIAL_CONDITION_LOAD_STEP){
		int NActives = MIN(INITIAL_CONDITION_LOAD_STEP,NProcs-i);
        for(int k=0;k<NActives;k++){
            if(MyID == i+k){
                fprintf(stderr,"Start read file. Node ID = %d\n",i+k);
                if(BinaryFlag){
                    ReadMusicLineBinary(fname,Numbers[0]+Numbers[1]+Numbers[2],mycount_comp,Eps,Params[7]);
                } else {
                    ReadMusicLine(fname,Numbers[0]+Numbers[1]+Numbers[2],mycount_comp,Eps,Params[7]);
                }
                fprintf(stderr,"Finish read file. Node ID = %d\n",i+k);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    CheckEpsValues();

    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.e0;

    Pall.Ns = 128;
    Pall.Npm = 8;

    Pall.hubble = Params[5];
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    Pall.OmegaM = Params[3];
    Pall.OmegaB = 0.0;
    Pall.OmegaL = Params[4];
    Pall.InitialRedshift = Pall.Redshift = Params[1];

    Pall.TEnd = CalcZtoT(0.e0);
    Pall.TCurrent = CalcZtoT(Pall.InitialRedshift);
    Pall.Redshift = Pall.InitialRedshift;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);
    UpdateAdaptiveSofteningFactor();

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = (Pall.TEnd-Pall.TCurrent)/5000.0; 
    Pall.OutPutInterval = (Pall.TEnd-Pall.TCurrent)/1000.0; 
    //Pall.OutPutInterval = 2.5*MEGAYEAR_CGS/Pall.UnitTime;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(!CheckDir("./data"))
            MakeDir("./data");
    }
    strcpy(Pall.ASCIIFileName,"./data/GF.ASCII");
    strcpy(Pall.BaseFileName,"./data/GF");
    strcpy(Pall.RestartFileName,"./data/GF.dump");

    return;
}

#endif // TASK_GALAXY_FORMATION  //}

struct StructHeader{
    float z_in;
    float H0;
    float Omega_m;
    float Omega_b;
    float Omega_l;
    float sigma_8;
    float n_s;
    float Lbox;
} Header;


void ReadMultiMassCosmologicalInitialCondition(char fbody[], char fboundary[]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = MPC_CGS;
    Pall.UnitTime = 10.0*GIGAYEAR_CGS;
    Pall.UnitMass = 1.e+11*MSUN_CGS;
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    // Set cosmological parameters!
    
    // read body! 
#if 0
    int count_body = 0;
    int count_boundary = 0;
    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,"SantaBarbara.Body.dat","r");
        while(fscanf(fp,"%*e %*e %*e %*e %*e %*e %*e") != EOF)
            count_body ++;
        fclose(fp);

        // read boundary! 
        FileOpen(fp,"SantaBarbara.Boundary.dat","r");
        while(fscanf(fp,"%*e %*e %*e %*e %*e %*e %*e") != EOF)
            count_boundary ++;
        fclose(fp);
    }
#endif

    int count_body = 0;
    int count_boundary = 0;
    for(int i=0;i<NProcs;i++){
        if(i == MyID){
            FILE *fp;
            FileOpen(fp,fbody,"rb");
            fread(&count_body,sizeof(int),1,fp);       
            fread(&count_boundary,sizeof(int),1,fp);       
            fread(&Header,sizeof(struct StructHeader),1,fp);       
            fclose(fp);
            fprintf(stderr,"[%02d] Body %d, Boundary %d\n ",MyID,count_body,count_boundary);
            fprintf(stderr,"z_in = %g\n",Header.z_in);
            fprintf(stderr,"H0 = %g\n",Header.H0);
            fprintf(stderr,"Omega_m = %g\n",Header.Omega_m);
            fprintf(stderr,"Omega_b = %g\n",Header.Omega_b);
            fprintf(stderr,"Omega_l = %g\n",Header.Omega_l);
            fprintf(stderr,"sigma_8 = %g\n",Header.sigma_8);
            fprintf(stderr,"n_s = %g\n",Header.n_s);
            fprintf(stderr,"Lbox = %g\n",Header.Lbox);
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /*
    if(count_body > 10000)
        count_body = 10000;
    if(count_boundary > 10000)
        count_boundary = 10000;
    */

    // MPI_Bcast(&count_body,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    // MPI_Bcast(&count_boundary,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);

    Pall.Ntotal_t = Pall.NDM_t = count_body + count_boundary;
    dlprintlmpi(Pall.Ntotal_t);
    dlprintlmpi(Pall.NDM_t);

    int mycount_body = 0;
    for(int k=0;k<count_body;k+=NProcs){
        int Element = MIN(count_body-k,NProcs);
        if(MyID < Element)
            mycount_body ++;
    }

    int mycount_boundary = 0;
    for(int k=0;k<count_boundary;k+=NProcs){
        int Element = MIN(count_boundary-k,NProcs);
        if(MyID < Element)
            mycount_boundary ++;
    }

    int mycount = mycount_body + mycount_boundary; 

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);
    //dprintlmpi(mycount);
    //dprintlmpi(AllocationSize);

    PbodySize = AllocationSize;
    PbodyElementsSize = AllocationSize;
    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(AllocationSize*sizeof(StructPbodyptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<AllocationSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Pbody[i] = PbodyElements+i;

    double EPS_CDM = 10*KPC_CGS/Pall.UnitLength;
    //double EPS_CDM = 5*KPC_CGS/Pall.UnitLength;
    for(int i=0;i<NProcs;i++){
        if(MyID == i){
            fprintf(stderr,"Read Data : %02d Node\n",i);
            int DummyInt;
            struct StructHeader DummyHeader;
            float PosVelMass[7];

            FILE *fp;
            FileOpen(fp,fbody,"rb");

            fread(&DummyInt,sizeof(int),1,fp);
            fread(&DummyInt,sizeof(int),1,fp);
            fread(&DummyHeader,sizeof(struct StructHeader),1,fp);
            int current = 0;
            for(int k=0;k<count_body;k++){
                fread(PosVelMass,sizeof(float),7,fp);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = PosVelMass[0];
                    Pbody[current]->Pos[1] = PosVelMass[1];
                    Pbody[current]->Pos[2] = PosVelMass[2];

                    Pbody[current]->Vel[0] = PosVelMass[3];
                    Pbody[current]->Vel[1] = PosVelMass[4];
                    Pbody[current]->Vel[2] = PosVelMass[5];

                    Pbody[current]->Mass = PosVelMass[6];

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    Pbody[current]->Type = TypeDM;

                    Pbody[current]->GlobalID = k;

                    current ++;
                }
            }
            fclose(fp);

            FileOpen(fp,fboundary,"rb");
            for(int k=0;k<count_boundary;k++){
                fread(PosVelMass,sizeof(float),7,fp);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = PosVelMass[0];
                    Pbody[current]->Pos[1] = PosVelMass[1];
                    Pbody[current]->Pos[2] = PosVelMass[2];

                    Pbody[current]->Vel[0] = PosVelMass[3];
                    Pbody[current]->Vel[1] = PosVelMass[4];
                    Pbody[current]->Vel[2] = PosVelMass[5];

                    Pbody[current]->Mass = PosVelMass[6];

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    Pbody[current]->Type = TypeDM;

                    //Pbody[current]->GlobalID = -1;
                    Pbody[current]->GlobalID = count_body+k;

                    current ++;
                }

            }
            fclose(fp);
            Pall.Ntotal = Pall.NDM = current;
            ReConnectPointers();
        }
    }
    UpdateTotalNumber();

    Pall.RunStatus = NewSimulation;

    Pall.hubble = 0.01*Header.H0;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    Pall.OmegaM = Header.Omega_m;
    Pall.OmegaL = Header.Omega_l;
    Pall.OmegaB = Header.Omega_b;
    Pall.InitialRedshift = Pall.Redshift = Header.z_in;
    fprintf(stderr,"h0 = %g, Om = %g, Ol = %g, Ob = %g, InitZ = %g\n",
            Pall.hubble,Pall.OmegaM,Pall.OmegaL,Pall.OmegaB,Pall.InitialRedshift);

    Pall.TEnd = CalcZtoT(0.e0);
    Pall.TCurrent = CalcZtoT(Pall.InitialRedshift);
    Pall.Redshift = Pall.InitialRedshift;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
            Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);
    UpdateAdaptiveSofteningFactor();

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/TheSmallBox.ASCII");
    strcpy(Pall.BaseFileName,"./data/TheSmallBox");
    strcpy(Pall.RestartFileName,"./data/TheSmallBox.dump");

    //exit(1);

    /*
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(stderr,"%g %g %g %g %g %g %g\n",
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],Pbody[i]->Mass);
    }
    exit(1);
    */
    OutPutAllParticlesInASCIIFormat();
    
    //InitLogFiles();

    return;
}

void ReadMultiMassNbodyInitialCondition(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = MPC_CGS;
    Pall.UnitTime = MPC_CGS/VELOCITY_KMS_CGS;
    Pall.UnitMass = 1.e+11*MSUN_CGS;
    Pall.TCMB = CMB_TEMPERATURE;

    eprintlmpi(Pall.UnitLength);
    eprintlmpi(Pall.UnitMass);
    eprintlmpi(Pall.UnitTime);

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 0.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    // Set cosmological parameters!
    
    // read body! 
    int count_body = 0;
    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,"Level0.dat","r");
        while(fscanf(fp,"%*e %*e %*e %*e %*e %*e %*e %*d %*d %*d %*u") != EOF)
            count_body ++;
        fclose(fp);
    }

    MPI_Bcast(&count_body,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);

    Pall.Ntotal_t = Pall.NDM_t = count_body;
    dlprintlmpi(Pall.Ntotal_t);
    dlprintlmpi(Pall.NDM_t);

    int mycount_body = 0;
    for(int k=0;k<count_body;k+=NProcs){
        int Element = MIN(count_body-k,NProcs);
        if(MyID < Element)
            mycount_body ++;
    }
    int mycount = mycount_body; 

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);
    dprintlmpi(mycount);
    dprintlmpi(AllocationSize);

    PbodySize = AllocationSize;
    PbodyElementsSize = AllocationSize;
    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(AllocationSize*sizeof(StructPbodyptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<AllocationSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Pbody[i] = PbodyElements+i;

    double EPS_CDM = 0.05;
    for(int i=0;i<NProcs;i++){
        if(MyID == i){

            double x,y,z,vx,vy,vz,mass;
            unsigned long int ID;

            FILE *fp;
            FileOpen(fp,"Level0.dat","r");
            int current = 0;
            for(int k=0;k<count_body;k++){
                fscanf(fp,"%le %le %le %le %le %le %le %*d %*d %*d %lu",
                        &x,&y,&z,&vx,&vy,&vz,&mass,&ID);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = x;
                    Pbody[current]->Pos[1] = y;
                    Pbody[current]->Pos[2] = z;

                    Pbody[current]->Vel[0] = vx;
                    Pbody[current]->Vel[1] = vy;
                    Pbody[current]->Vel[2] = vz;

                    Pbody[current]->Mass = mass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    Pbody[current]->Type = TypeDM;

                    Pbody[current]->GlobalID = ID;

                    current ++;
                }

            }
            fclose(fp);


            Pall.Ntotal = Pall.NDM = current;

            ReConnectPointers();
            UpdateTotalNumber();
        }
    }

    Pall.RunStatus = NewSimulation;

    FILE *fp;
    // Read cosmological parameters.
    FileOpen(fp,"CosmologicalParameters.dat","r");
    fscanf(fp,"%le %le %le",&Pall.OmegaB,&Pall.OmegaM,&Pall.OmegaL);
    fscanf(fp,"%le %*g %*g",&Pall.hubble);
    fclose(fp);
    Pall.OmegaM = 1.e0 - Pall.OmegaL;
    Pall.hubble *= 0.01;
    fprintf(stderr,"Cosmological Parameters : (OmegaB, OmegaCDM, OmegaL, Hubble) = (%g, %g, %g, %g)\n",
            Pall.OmegaB,Pall.OmegaM,Pall.OmegaL,Pall.hubble);

    // Read redshift etc.
    FileOpen(fp,"grid.dx.astart","r");
    fscanf(fp,"%*d %*d %*d %*g %le",&Pall.InitialRedshift);
    fclose(fp);
    Pall.Redshift = Pall.InitialRedshift = 1.0/Pall.InitialRedshift - 1.0;
    fprintf(stderr,"Scale Factor = %g, Initial Redshift = %g\n",1.0/(1+Pall.InitialRedshift),Pall.InitialRedshift);


    //Pall.hubble = 0.5;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    //Pall.OmegaM = 1.0;
    //Pall.OmegaL = 0.0;
    //Pall.OmegaB = 0.0;
    //Pall.InitialRedshift = Pall.Redshift = 63.0;

    Pall.TEnd = CalcZtoT(0.e0);
    Pall.TCurrent = CalcZtoT(Pall.InitialRedshift);
    Pall.Redshift = Pall.InitialRedshift;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
            Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    fprintf(stderr,"Tstart(in Gyr) = %g, Tend(in Gyr) = %g\n",
            Pall.TCurrent*Pall.UnitTime/GIGAYEAR_CGS,Pall.TEnd*Pall.UnitTime/GIGAYEAR_CGS);
    UpdateAdaptiveSofteningFactor();

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",Pall.AdaptiveSofteningFactor);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1*GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/Nbody.ASCII");
    strcpy(Pall.BaseFileName,"./data/Nbody");
    strcpy(Pall.RestartFileName,"./data/Nbody.dump");


    InitLogFiles();

    return;
}

#ifdef TASK_NFW
double M200_NFW;
double Tdyn_NFW;
double Cnfw_NFW;
double Spin_NFW;
double fbaryon_NFW;

//inline double fx(const double x) __attribute__((always_inline));
//inline double fx(const double x){
static inline double __attribute__((always_inline)) fx(const double x){
    return (log(1.0+x) - x/(1.e0+x));
}

//void ExternalPotentialForNFW(const double M200){
void ExternalPotentialForNFW(void){

    //double M200 = pow(10.0,(double)(M200_Log10M-10));
    //dprintlmpi(M200_Log10M);
    //gprintlmpi(M200);
    //exit(0);
    //double Cnfw = 10.e0;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double delta_c = 200 * CUBE(Cnfw_NFW)/(3*fx(Cnfw_NFW));
    double rho0 = RhoCrit * delta_c;
    double R200 = pow(M200_NFW*3/(200*RhoCrit*4*M_PI), 1.e0/3.e0);
    double rs = R200/Cnfw_NFW;

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            double r = sqrt(NORM2(Pbody[i]->PosP)+SQ(0.1*rs));
            double x = r/rs;
            double Mass = 4*M_PI*rho0*CUBE(rs)*fx(x);
            double gravfact = Pall.GravConst*Mass/CUBE(r);

            Pbody[i]->Acc[0] -= gravfact*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] -= gravfact*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] -= gravfact*Pbody[i]->PosP[2];
        }
    }

    return ;
}

void ReadNFWInitialCondition(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MPIGetMyID());

    double M200,Tdyn;
    double Cnfw,fbaryon,Spin;

    //double Tdyn,epsilon,M200;
    double epsilon;
    for(int k=0;k<NProcs;k++){
        if(k == MyID){
            FILE *fp;
            char fname[MaxCharactersInLine];
            sprintf(fname, "ASURA_nfw.dat");
            fprintf(stderr,"Read Initial Constion[%d] from %s\n",MyID,fname);
            FileOpen(fp,fname,"rb");

            int ReadNumber;
            fread(&ReadNumber,sizeof(int),1,fp);
            dprintlmpi(ReadNumber);
            Pall.Ntotal_t = ReadNumber;
            Pall.Nhydro_t = ReadNumber;

            double Units;
            fread(&Units,sizeof(double),1,fp);
            Pall.UnitLength = Units;
            fread(&Units,sizeof(double),1,fp);
            Pall.UnitMass = Units;
            fread(&Units,sizeof(double),1,fp);
            Pall.UnitTime = Units;

            fread(&Tdyn,sizeof(double),1,fp);
            fread(&M200,sizeof(double),1,fp);
            fread(&epsilon,sizeof(double),1,fp);
            M200_NFW = M200;
            Tdyn_NFW = Tdyn;

            fread(&Cnfw,sizeof(double),1,fp);
            fread(&Spin,sizeof(double),1,fp);
            fread(&fbaryon,sizeof(double),1,fp);
            Cnfw_NFW = Cnfw;
            Spin_NFW = Spin;
            fbaryon_NFW = fbaryon;

            fclose(fp);

            if(k==0){
                fprintf(stderr,"Number = %d, L = %g, M = %g T = %g, Tdyn = %g",
                    ReadNumber,Pall.UnitLength,Pall.UnitMass,Pall.UnitTime,Tdyn);
                fprintf(stderr,"M200 = %g, eps = %g, Cnfw = %g, Spin = %g, fbaryon = %g\n",
                    M200,epsilon,Cnfw,Spin,fbaryon);

                fprintf(stderr,"%d %g kpc, %g Msun, %g Year, Tdyn = %g Year, eps = %g kpc\n",
                    ReadNumber,Pall.UnitLength/KPC_CGS,Pall.UnitMass/MSUN_CGS,
                        Pall.UnitTime/YEAR_CGS,
                        Pall.UnitTime*Tdyn/YEAR_CGS,Pall.UnitLength*epsilon/KPC_CGS);
            }

            for(int i=0;i<Pall.Ntotal_t;i++){
                if(i%NProcs == MyID){
                    Pall.Nhydro ++;
                    Pall.Ntotal ++;
                }
            }
        }
        // MPI_Barrier(MPI_COMM_WORLD);
    }


    //Pall.UnitLength = MPC_CGS;
    //Pall.UnitTime = 10.0*GIGAYEAR_CGS;
    //Pall.UnitMass = 1.e+11*MSUN_CGS;
    //Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 10.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}


    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    PbodySize = Pall.Ntotal;
    PhydroSize = Pall.Nhydro;
    PbodyElementsSize = Pall.Ntotal;
    PhydroElementsSize = Pall.Nhydro;

    dlprintlmpi(Pall.NDM);
    dlprintlmpi(Pall.Ntotal);
    dlprintlmpi(Pall.Nhydro);


    GenerateStructPbody(Pall.Ntotal);
    GenerateStructPhydro(Pall.Nhydro);

    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    struct StructInput{
        float Pos[3];
        float Vel[3];
        float Mass;
        float EgySpec;
        int Id;
    } InputBasket;

    //double Eps = 10*PC_CGS/Pall.UnitLength;
    double Eps = 25*PC_CGS/Pall.UnitLength;

    for(int k=0;k<NProcs;k++){
        if(k == MyID){
            FILE *fp;
            char fname[MaxCharactersInLine];
            sprintf(fname, "ASURA_nfw.dat");
            FileOpen(fp,fname,"rb");
            fseek(fp,sizeof(int)+9*sizeof(double),SEEK_SET);

            // dprintlmpi(Pall.Ntotal_t);

            int mycount = 0;
            for(int i=0;i<Pall.Ntotal_t;i++){
                fread(&InputBasket,sizeof(struct StructInput),1,fp);
                if(i%NProcs == MyID){
                    Pbody[mycount]->Active = ON;
                    Pbody[mycount]->Use = ON;
                    Pbody[mycount]->Type = TypeHydro;
                    Pbody[mycount]->GlobalID = (double)InputBasket.Id;

                    Pbody[mycount]->Pos[0] = (double)InputBasket.Pos[0];
                    Pbody[mycount]->Pos[1] = (double)InputBasket.Pos[1];
                    Pbody[mycount]->Pos[2] = (double)InputBasket.Pos[2];
                    Pbody[mycount]->PosP[0] = Pbody[mycount]->Pos[0];
                    Pbody[mycount]->PosP[1] = Pbody[mycount]->Pos[1];
                    Pbody[mycount]->PosP[2] = Pbody[mycount]->Pos[2];

                    Pbody[mycount]->Vel[0] = (double)InputBasket.Vel[0];
                    Pbody[mycount]->Vel[1] = (double)InputBasket.Vel[1];
                    Pbody[mycount]->Vel[2] = (double)InputBasket.Vel[2];

                    Pbody[mycount]->Mass = (double)InputBasket.Mass;
                    //Pbody[mycount]->Eps = epsilon;
                    Pbody[mycount]->Eps = Eps;

                    PbodyHydro(mycount)->Use = ON;
                    PbodyHydro(mycount)->Kernel = 2.0*epsilon;
                    PbodyHydro(mycount)->U = (double)InputBasket.EgySpec;

#if (UseSFModelSpawn) 
                    PbodyHydro(mycount)->SpawnMass = Pbody[mycount]->Mass/(double)MaxSpawnTimes;
#endif
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
                    PbodyHydro(mycount)->G0 = 1;
                    PbodyHydro(mycount)->fH2 = 0.1;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS


#ifdef USE_CELIB //{
                    CELibSetPrimordialMetallicity(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements);
                    //CELibSetSolarMetallicity(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements);
                    // CELibSetMetallicityWithSolarAbundancePattern(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements,
                            // 0.1*CELibGetMetalFractionForSolarChemicalComposision());
                    double MassLightElements = 
                        Phydro[mycount]->Elements[CELibYield_H]+Phydro[mycount]->Elements[CELibYield_He];
                    double Z = (Pbody[mycount]->Mass-MassLightElements)/Pbody[mycount]->Mass;
                    PbodyHydro(mycount)->Z = PbodyHydro(mycount)->ZII = PbodyHydro(mycount)->ZIa = Z;
#else 
                    PbodyHydro(mycount)->Z   = 0.02;
                    PbodyHydro(mycount)->ZII = 0.02;
                    PbodyHydro(mycount)->ZIa = 0.00;
#endif // USE_CELIB //}

                    mycount ++;
                }
            }
            fclose(fp);
        }
        // MPI_Barrier(MPI_COMM_WORLD);
    }


    ActivateAllparticles();
#if 0
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }
#endif

#if 1
    FILE *fp;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"Init.%04d.%04d",NProcs,MyID);
    FileOpen(fp ,fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pall.ConvertUtoT*PbodyHydro(i)->U,Pbody[i]->Mass);
    }
    fclose(fp);
    fflush(NULL);

    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat ./Init.????.???? | sort -n > ./Init.dat");
        fflush(NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("rm -rf ./Init.????.????");
        fflush(NULL);
    }
    //exit(0);
#endif

    // ReConnectPointers();
    // UpdateTotalNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.e0;

    Pall.Ns = 128;
    Pall.Npm = 8;

    Pall.hubble = 0.01*70;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    //Pall.OmegaM = 1.0;
    //Pall.OmegaL = 0.0; 

    Pall.TEnd = 5*GIGAYEAR_CGS/Pall.UnitTime;
    // Pall.TEnd = 5*MEGAYEAR_CGS/Pall.UnitTime;
    //Pall.TEnd = 1*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = 0.1*(Tdyn); 
    Pall.OutPutInterval = Pall.TEnd/1000.0; 
    //Pall.OutPutInterval = 10*Pall.TEnd; 

    if(MPIGetMyID() == MPI_ROOT_RANK)
        MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/NFW.ASCII");
    strcpy(Pall.BaseFileName,"./data/NFW");
    strcpy(Pall.RestartFileName,"./data/NFW.dump");

    return;
}


void ReadNFWTdyn(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    double M200,Tdyn;
    double epsilon;
    double Cnfw,fbaryon,Spin;

    for(int k=0;k<NProcs;k++){
        if(k == MyID){
            FILE *fp;
            char fname[MaxCharactersInLine];
            sprintf(fname, "ASURA_nfw.dat");
            fprintf(stderr,"Read Initial Constion[%d] from %s\n",MyID,fname);
            FileOpen(fp,fname,"rb");

            int ReadNumber;
            fread(&ReadNumber,sizeof(int),1,fp);
            double Units;
            fread(&Units,sizeof(double),1,fp);
            fread(&Units,sizeof(double),1,fp);
            fread(&Units,sizeof(double),1,fp);

            fread(&Tdyn,sizeof(double),1,fp);
            fread(&M200,sizeof(double),1,fp);
            fread(&epsilon,sizeof(double),1,fp);
            M200_NFW = M200;
            Tdyn_NFW = Tdyn;

            fread(&Cnfw,sizeof(double),1,fp);
            fread(&Spin,sizeof(double),1,fp);
            fread(&fbaryon,sizeof(double),1,fp);
            Cnfw_NFW = Cnfw;
            Spin_NFW = Spin;
            fbaryon_NFW = fbaryon;

            fclose(fp);

            if(k==0){
                fprintf(stderr,"Number = %d, L = %g, M = %g T = %g, Tdyn = %g",
                    ReadNumber,Pall.UnitLength,Pall.UnitMass,Pall.UnitTime,Tdyn);
                fprintf(stderr,"M200 = %g, eps = %g, Cnfw = %g, Spin = %g, fbaryon = %g\n",
                    M200,epsilon,Cnfw,Spin,fbaryon);

                fprintf(stderr,"%d %g kpc, %g Msun, %g Year, Tdyn = %g Year, eps = %g kpc\n",
                    ReadNumber,Pall.UnitLength/KPC_CGS,Pall.UnitMass/MSUN_CGS,
                        Pall.UnitTime/YEAR_CGS,Tdyn*Pall.UnitTime/YEAR_CGS,epsilon*Pall.UnitLength/KPC_CGS);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    return; 
}
#endif

static double ReturnDistanceForExponent(const double F, const double Rd){

#define Rstart (0)
#define Rend   (5*Rd)
//#define Rend   (0.5*Rd)

    double Diff;

    double rstart = Rstart;
    double rend = Rend;
    double fstart = 1-exp(-rstart/Rd)*(1+rstart/Rd) - F;
    double fend = 1-exp(-rend/Rd)*(1+rend/Rd) - F;

    //gprintlmpi(fstart);
    //gprintlmpi(fend);

    double r = 0.5*(rend+rstart);
    double f = 1-exp(-r/Rd)*(1+r/Rd) - F;

    //gprintlmpi(r);
    //gprintlmpi(f);

    //if(F > fend)
        //return Rend;

    int Iteration = 0;
    do{
        f = 1-exp(-r/Rd)*(1+r/Rd) - F;

        if(f*fstart < 0.0){
            fend = f;
            rend = r;
            r = 0.5*(rstart+rend);
            Diff = fabs(rend-rstart);
        } else {
            fstart = f;
            rstart = r;
            r = 0.5*(rstart+rend);
            Diff = fabs(rend-rstart);
        }
        //gprintlmpi(Diff);
        if(Iteration > MaxIteration)
            exit(1);

        Iteration ++;
    } while ( Diff > 1.e-6*Rd );

#undef Rstart 
#undef Rend  
    return r;
}

static double ReturnDistanceForInvR(const double F, const double Rd){

    double Rdstart = (0.0);
    double Rdend   = (10.0*Rd);

    double Diff;

    double rstart = Rdstart;
    double rend = Rdend;
    double fstart = rstart/Rd - F;
    double fend = rend/Rd - F;

    //gprintlmpi(fstart);
    //gprintlmpi(fend);

    double r = 0.5*(rend+rstart);
    double f = r/Rd - F;

    //gprintlmpi(r);
    //gprintlmpi(f);

    //if(F > fend)
        //return Rend;

    int Iteration = 0;
    do{
        f = r/Rd - F;

        if(f*fstart < 0.0){
            fend = f;
            rend = r;
            r = 0.5*(rstart+rend);
            Diff = fabs(rend-rstart);
        } else {
            fstart = f;
            rstart = r;
            r = 0.5*(rstart+rend);
            Diff = fabs(rend-rstart);
        }
        //gprintlmpi(Diff);
        if(Iteration > MaxIteration)
            exit(1);

        Iteration ++;
    } while ( Diff > 1.e-6*Rd );

    return r;
}

#ifdef TASK_MERGER
void ReadGalactICS(void){

#define GalactICS_UnitMass (5.e+10)     //5x10^10 Msun
#define GalactICS_UnitLength (4.5)      //4.5 kpc
#define GalactICS_UnitVelocity (220.0)  //220 km/s

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = KPC_CGS;
    Pall.UnitMass = MSUN_CGS;
    Pall.UnitTime = KPC_CGS/(VELOCITY_KMS_CGS);
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    // Set cosmological parameters!
    
    // read body! 
    int count_disk,count_bulge,count_halo;
    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,"./disk","r");
        fscanf(fp,"%d %*e",&count_disk);
        fclose(fp);

        FileOpen(fp,"./bulge","r");
        fscanf(fp,"%d %*e",&count_bulge);
        fclose(fp);

        FileOpen(fp,"./halo","r");
        fscanf(fp,"%d %*e",&count_halo);
        fclose(fp);
    }

    MPI_Bcast(&count_disk,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&count_bulge,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&count_halo,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);

    Pall.Ntotal_t = Pall.NDM_t = count_disk + count_bulge + count_halo;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Ntotal_t = %d, NDM_t = %d\n",Pall.Ntotal_t,Pall.NDM_t);

    int mycount_disk = 0;
    for(int k=0;k<count_disk;k+=NProcs){
        int Element = MIN(count_disk-k,NProcs);
        if(MyID < Element)
            mycount_disk ++;
    }
    int mycount_bulge = 0;
    for(int k=0;k<count_bulge;k+=NProcs){
        int Element = MIN(count_bulge-k,NProcs);
        if(MyID < Element)
            mycount_bulge ++;
    }
    int mycount_halo = 0;
    for(int k=0;k<count_halo;k+=NProcs){
        int Element = MIN(count_halo-k,NProcs);
        if(MyID < Element)
            mycount_halo ++;
    }
    int mycount = mycount_disk + mycount_bulge + mycount_halo; 
    Pall.Ntotal = Pall.NDM = mycount;
    fprintf(stderr,"mycount_disk = %d, mycount_bulge = %d, mycount_halo = %d\n",mycount_disk,mycount_bulge,mycount_halo);
    fprintf(stderr,"Ntotal = %d, NDM = %d\n",Pall.Ntotal,Pall.NDM);

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);

    PbodySize = AllocationSize;
    PbodyElementsSize = AllocationSize;
    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(AllocationSize*sizeof(StructPbodyptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<AllocationSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);
    PbodyElements[AllocationSize-1].Next = NULL;

    for(int i=0;i<AllocationSize;i++)
        Pbody[i] = PbodyElements+i;

    //double EPS_CDM = 200*PC_CGS/Pall.UnitLength;
    double EPS_CDM = 20*PC_CGS/Pall.UnitLength;
    //double EPS_CDM = 10*PC_CGS/Pall.UnitLength;
    if(MyID == MPI_ROOT_RANK)
        fprintf(stderr,"The softening length = %g pc\n",EPS_CDM*Pall.UnitLength/PC_CGS);

    for(int i=0;i<NProcs;i++){
        if(MyID == i){

            double x,y,z,vx,vy,vz,mass;

            FILE *fp;
            FileOpen(fp,"./disk","r");
            int current = 0;
            fscanf(fp,"%*d %*e");
            for(int k=0;k<count_disk;k++){
                fscanf(fp,"%le %le %le %le %le %le %le",&mass,&x,&y,&z,&vx,&vy,&vz);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = x*GalactICS_UnitLength;
                    Pbody[current]->Pos[1] = y*GalactICS_UnitLength;
                    Pbody[current]->Pos[2] = z*GalactICS_UnitLength;

                    Pbody[current]->Vel[0] = vx*GalactICS_UnitVelocity;
                    Pbody[current]->Vel[1] = vy*GalactICS_UnitVelocity;
                    Pbody[current]->Vel[2] = vz*GalactICS_UnitVelocity;

                    Pbody[current]->Mass = mass*GalactICS_UnitMass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    Pbody[current]->Type = TypeDM;

                    Pbody[current]->GlobalID = k;

                    current ++;
                }

            }
            fclose(fp);


            FileOpen(fp,"./bulge","r");
            fscanf(fp,"%*d %*e");
            for(int k=0;k<count_bulge;k++){
                fscanf(fp,"%le %le %le %le %le %le %le",&mass,&x,&y,&z,&vx,&vy,&vz);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = x*GalactICS_UnitLength;
                    Pbody[current]->Pos[1] = y*GalactICS_UnitLength;
                    Pbody[current]->Pos[2] = z*GalactICS_UnitLength;

                    Pbody[current]->Vel[0] = vx*GalactICS_UnitVelocity;
                    Pbody[current]->Vel[1] = vy*GalactICS_UnitVelocity;
                    Pbody[current]->Vel[2] = vz*GalactICS_UnitVelocity;

                    Pbody[current]->Mass = mass*GalactICS_UnitMass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    Pbody[current]->Type = TypeDM;

                    Pbody[current]->GlobalID = count_disk + k;

                    current ++;
                }

            }
            fclose(fp);

            FileOpen(fp,"./halo","r");
            fscanf(fp,"%*d %*e");
            for(int k=0;k<count_halo;k++){
                fscanf(fp,"%le %le %le %le %le %le %le",&mass,&x,&y,&z,&vx,&vy,&vz);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = x*GalactICS_UnitLength;
                    Pbody[current]->Pos[1] = y*GalactICS_UnitLength;
                    Pbody[current]->Pos[2] = z*GalactICS_UnitLength;

                    Pbody[current]->Vel[0] = vx*GalactICS_UnitVelocity;
                    Pbody[current]->Vel[1] = vy*GalactICS_UnitVelocity;
                    Pbody[current]->Vel[2] = vz*GalactICS_UnitVelocity;

                    Pbody[current]->Mass = mass*GalactICS_UnitMass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    Pbody[current]->Type = TypeDM;

                    Pbody[current]->GlobalID = count_disk + count_bulge + k;

                    current ++;
                }

            }
            fclose(fp);
            dprintlmpi(current);
            dprintlmpi(Pall.Ntotal);
            assert(current == Pall.Ntotal);
            Pall.Ntotal = Pall.NDM = current;

        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    ReConnectPointers();
    UpdateTotalNumber();

    Pall.RunStatus = NewSimulation;

    // Pall.hubble = 0.5;
    // Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    // Pall.OmegaM = 1.0;
    // Pall.OmegaL = 0.0;
    // Pall.OmegaB = 0.0;
    Pall.InitialRedshift = Pall.Redshift = 0;

    Pall.TEnd = 2.0*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g\n",Pall.TCurrent,Pall.TEnd);

    Pall.AdaptiveSofteningFactor = 1.e0;
    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1*GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/ID.ASCII");
    strcpy(Pall.BaseFileName,"./data/ID");
    strcpy(Pall.RestartFileName,"./data/ID.dump");

    //InitLogFiles();

    OutPutDarkMatter("init.dat");

    return;
}

struct StructReadGalactICSBody{
    float Mass;
    float Pos[3];
    float Vel[3];
};

void ReadBinaryGalactICS(char fname[], const double GasFlactionDisk, const double GasFlactionHalo){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = KPC_CGS;
    Pall.UnitMass = MSUN_CGS;
    Pall.UnitTime = KPC_CGS/(VELOCITY_KMS_CGS);
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    // Set cosmological parameters!
    
    // read body! 
    int count_disk,count_bulge,count_halo;
    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,fname,"rb");
        fread(&count_disk,sizeof(int),1,fp);
        fread(&count_bulge,sizeof(int),1,fp);
        fread(&count_halo,sizeof(int),1,fp);
        fclose(fp);
    }
    MPI_Bcast(&count_disk,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&count_bulge,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&count_halo,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);

    // hydro
    int count_disk_gas = count_disk*GasFlactionDisk; 
    int count_halo_gas = count_halo*GasFlactionHalo; 
    //MPI_Bcast(&count_disk_gas,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    //MPI_Bcast(&count_disk_halo,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    fprintf(stderr,"[%02d] Ndisk,Nbulge,Nhalo = %d %d %d\n",MyID,count_disk,count_bulge,count_halo);
    fprintf(stderr,"[%02d] Ndisk_gas,Nhalo_gas = %d %d\n",MyID,count_disk_gas,count_halo_gas);

    Pall.Ntotal_t = Pall.NDM_t = count_disk + count_bulge + count_halo;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Ntotal_t = %d, NDM_t = %d\n",Pall.Ntotal_t,Pall.NDM_t);

    int mycount_disk = 0;
    for(int k=0;k<count_disk;k+=NProcs){
        int Element = MIN(count_disk-k,NProcs);
        if(MyID < Element)
            mycount_disk ++;
    }
    int mycount_bulge = 0;
    for(int k=0;k<count_bulge;k+=NProcs){
        int Element = MIN(count_bulge-k,NProcs);
        if(MyID < Element)
            mycount_bulge ++;
    }
    int mycount_halo = 0;
    for(int k=0;k<count_halo;k+=NProcs){
        int Element = MIN(count_halo-k,NProcs);
        if(MyID < Element)
            mycount_halo ++;
    }
    int mycount = mycount_disk + mycount_bulge + mycount_halo; 
    Pall.Ntotal = Pall.NDM = mycount;
    fprintf(stderr,"mycount_disk = %d, mycount_bulge = %d, mycount_halo = %d\n",mycount_disk,mycount_bulge,mycount_halo);
    fprintf(stderr,"Ntotal = %d, NDM = %d\n",Pall.Ntotal,Pall.NDM);

    // hydro
    int mycount_disk_gas = mycount_disk*GasFlactionDisk; 
    int mycount_halo_gas = mycount_halo*GasFlactionHalo; 
    Pall.Nhydro = mycount_disk_gas+mycount_halo_gas;
    int LocalNhydro,GlobalNhydro;
    LocalNhydro = Pall.Nhydro;
    MPI_Allreduce(&LocalNhydro,&GlobalNhydro,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    Pall.Nhydro_t = GlobalNhydro;

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);
    // Connect Body and Hydro locally.
    for(int k=0;k<Pall.Nhydro;k++){
        Pbody[k] = PbodyElements+k;
        Phydro[k] = PhydroElements+k;
        Pbody[k]->Baryon = (void *)(Phydro[k]);
        Phydro[k]->Body = Pbody[k];
    }

    //PbodySize = AllocationSize;
    //PbodyElementsSize = AllocationSize;
    //PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    //Pbody = malloc(AllocationSize*sizeof(StructPbodyptr));
    //memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));
//
    //for(int i=0;i<AllocationSize-1;i++)
        //PbodyElements[i].Next = &(PbodyElements[i+1]);
    //PbodyElements[AllocationSize-1].Next = NULL;
//
    //for(int i=0;i<AllocationSize;i++)
        //Pbody[i] = PbodyElements+i;

    //double EPS_CDM = 200*PC_CGS/Pall.UnitLength;
    double EPS_CDM = 20*PC_CGS/Pall.UnitLength;
    //double EPS_CDM = 10*PC_CGS/Pall.UnitLength;
    if(MyID == MPI_ROOT_RANK)
        fprintf(stderr,"The softening length = %g pc\n",EPS_CDM*Pall.UnitLength/PC_CGS);
    // O.K.


    for(int i=0;i<NProcs;i++){
        if(MyID == i){

            struct StructReadGalactICSBody Body;

            FILE *fp;
            FileOpen(fp,fname,"rb");
            fseek(fp,sizeof(int)*3,SEEK_SET);

            int current_hydro = 0;
            int current_nbody = 0;

            int current = 0;
            for(int k=0;k<count_disk;k++){
                fread(&Body,sizeof(struct StructReadGalactICSBody),1,fp);
                if(k%NProcs == MyID){
                    if(current_hydro < mycount_disk_gas){
                        int index = current_hydro; 
                        Pbody[index]->Pos[0] = Body.Pos[0];
                        Pbody[index]->Pos[1] = Body.Pos[1];
                        Pbody[index]->Pos[2] = Body.Pos[2];

                        Pbody[index]->Vel[0] = Body.Vel[0];
                        Pbody[index]->Vel[1] = Body.Vel[1];
                        Pbody[index]->Vel[2] = Body.Vel[2];

                        Pbody[index]->Mass = Body.Mass;

                        Pbody[index]->Active = ON;
                        Pbody[index]->Use = ON;
                        Pbody[index]->Eps = EPS_CDM;
                        Pbody[index]->Type = TypeDM;

                        Pbody[index]->GlobalID = k;

                        current_hydro ++;
                    } else {
                        int index = current_nbody + Pall.Nhydro; 
                        Pbody[index]->Pos[0] = Body.Pos[0];
                        Pbody[index]->Pos[1] = Body.Pos[1];
                        Pbody[index]->Pos[2] = Body.Pos[2];

                        Pbody[index]->Vel[0] = Body.Vel[0];
                        Pbody[index]->Vel[1] = Body.Vel[1];
                        Pbody[index]->Vel[2] = Body.Vel[2];

                        Pbody[index]->Mass = Body.Mass;

                        Pbody[index]->Active = ON;
                        Pbody[index]->Use = ON;
                        Pbody[index]->Eps = EPS_CDM;
                        Pbody[index]->Type = TypeDM;

                        Pbody[index]->GlobalID = k;
                        current_nbody ++;
                    }
                    current ++;
                }
            }


            for(int k=0;k<count_bulge;k++){
                fread(&Body,sizeof(struct StructReadGalactICSBody),1,fp);
                if(k%NProcs == MyID){
                    int index = current_nbody + Pall.Nhydro; 
                    Pbody[index]->Pos[0] = Body.Pos[0];
                    Pbody[index]->Pos[1] = Body.Pos[1];
                    Pbody[index]->Pos[2] = Body.Pos[2];

                    Pbody[index]->Vel[0] = Body.Vel[0];
                    Pbody[index]->Vel[1] = Body.Vel[1];
                    Pbody[index]->Vel[2] = Body.Vel[2];

                    Pbody[index]->Mass = Body.Mass;

                    Pbody[index]->Active = ON;
                    Pbody[index]->Use = ON;
                    Pbody[index]->Eps = EPS_CDM;
                    Pbody[index]->Type = TypeDM;

                    Pbody[index]->GlobalID = count_disk + k;

                    current_nbody ++;
                    current ++;
                }

            }

            for(int k=0;k<count_halo;k++){
                fread(&Body,sizeof(struct StructReadGalactICSBody),1,fp);
                if(k%NProcs == MyID){
                    if(current_hydro < mycount_disk_gas+mycount_halo_gas){
                        int index = current_hydro; 
                        Pbody[index]->Pos[0] = Body.Pos[0];
                        Pbody[index]->Pos[1] = Body.Pos[1];
                        Pbody[index]->Pos[2] = Body.Pos[2];

                        Pbody[index]->Vel[0] = Body.Vel[0];
                        Pbody[index]->Vel[1] = Body.Vel[1];
                        Pbody[index]->Vel[2] = Body.Vel[2];

                        Pbody[index]->Mass = Body.Mass;

                        Pbody[index]->Active = ON;
                        Pbody[index]->Use = ON;
                        Pbody[index]->Eps = EPS_CDM;
                        Pbody[index]->Type = TypeDM;

                        Pbody[index]->GlobalID = count_disk + count_bulge + k;
                        current_hydro ++;
                    } else {
                        int index = current_nbody + Pall.Nhydro; 
                        Pbody[index]->Pos[0] = Body.Pos[0];
                        Pbody[index]->Pos[1] = Body.Pos[1];
                        Pbody[index]->Pos[2] = Body.Pos[2];

                        Pbody[index]->Vel[0] = Body.Vel[0];
                        Pbody[index]->Vel[1] = Body.Vel[1];
                        Pbody[index]->Vel[2] = Body.Vel[2];

                        Pbody[index]->Mass = Body.Mass;

                        Pbody[index]->Active = ON;
                        Pbody[index]->Use = ON;
                        Pbody[index]->Eps = EPS_CDM;
                        Pbody[index]->Type = TypeDM;

                        Pbody[index]->GlobalID = count_disk + count_bulge + k;

                        current_nbody ++;
                    }
                    current ++;
                }
            }
            fclose(fp);
            dprintlmpi(current);
            dprintlmpi(current_hydro);
            dprintlmpi(current_nbody);
            dprintlmpi(Pall.Ntotal);
            dprintlmpi(Pall.Nhydro);
            dprintlmpi(Pall.Ntotal_t);
            dprintlmpi(Pall.Nhydro_t);
            assert(current == Pall.Ntotal);
            assert(current_hydro == Pall.Nhydro);
            Pall.Ntotal = Pall.NDM = current;

        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    // mycount_disk_gas -> Gas; Change type and add extra information.
    // mycount_halo_gas -> Gas; Cahnge type and add extra information.

    for(int i=0;i<mycount_disk_gas+mycount_halo_gas;i++){
        PhydroBody(i)->Type = TypeHydro;

        Phydro[i]->Use = ON;
        Phydro[i]->Active = ON;
        if(i<mycount_disk_gas){
            Phydro[i]->U = Pall.ConvertTtoU*1.e+4;
        } else {
            Phydro[i]->U = Pall.ConvertTtoU*1.e+5;
        }
        Phydro[i]->Kernel = 10*PhydroBody(i)->Eps;
#if (UseSFModelSpawn) 
        Phydro[i]->SpawnMass = PhydroBody(i)->Mass/(double)MaxSpawnTimes;
#endif
        Phydro[i]->Z = Phydro[i]->ZII = Phydro[i]->ZIa = 0.01; 
        // We here assume the half solar metallicity for the ISM at initially.

        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
    }

    fprintf(stderr,"Thalo = %g K.\n",
            (0.59*PROTON_MASS_CGS/(3.0*BOLTZMANN_CONSTANT_CGS))*(GRAVITY_CONSTANT_CGS*1.2e+11*MSUN_CGS/(100*KPC_CGS)));

    ReConnectPointers();
    UpdateTotalNumber();


    //MPI_Barrier(MPI_COMM_WORLD);
    //exit(1);

    Pall.RunStatus = NewSimulation;
    // hydro parameters
    Pall.Ns = 32;
    Pall.Npm = 2;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    // Pall.hubble = 0.5;
    // Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    // Pall.OmegaM = 1.0;
    // Pall.OmegaL = 0.0;
    // Pall.OmegaB = 0.0;
    Pall.InitialRedshift = Pall.Redshift = 0;

    Pall.TEnd = 3.0*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g\n",Pall.TCurrent,Pall.TEnd);

    Pall.AdaptiveSofteningFactor = 1.e0;
    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1*GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/ID.ASCII");
    strcpy(Pall.BaseFileName,"./data/ID");
    strcpy(Pall.RestartFileName,"./data/ID.dump");


    return;
}

struct StructReadGalactICSBodyFlag{
    bool GasFlag;
    float Mass;
    float Pos[3];
    float Vel[3];
};

void ReadBinaryGalactICSFlag(char fname[]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = KPC_CGS;
    Pall.UnitMass = MSUN_CGS;
    Pall.UnitTime = KPC_CGS/(VELOCITY_KMS_CGS);
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    // Set cosmological parameters!
    
    // read body! 
    int count_disk,count_bulge,count_halo;
    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,fname,"rb");
        fread(&count_disk,sizeof(int),1,fp);
        fread(&count_bulge,sizeof(int),1,fp);
        fread(&count_halo,sizeof(int),1,fp);
        fclose(fp);
    }
    MPI_Bcast(&count_disk,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&count_bulge,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&count_halo,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);

    fprintf(stderr,"[%02d] Ndisk,Nbulge,Nhalo = %d %d %d\n",MyID,count_disk,count_bulge,count_halo);

    Pall.Ntotal_t = Pall.NDM_t = count_disk + count_bulge + count_halo;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Ntotal_t = %d, NDM_t = %d\n",Pall.Ntotal_t,Pall.NDM_t);

    int mycount_disk = 0;
    for(int k=0;k<count_disk;k+=NProcs){
        int Element = MIN(count_disk-k,NProcs);
        if(MyID < Element)
            mycount_disk ++;
    }
    int mycount_bulge = 0;
    for(int k=0;k<count_bulge;k+=NProcs){
        int Element = MIN(count_bulge-k,NProcs);
        if(MyID < Element)
            mycount_bulge ++;
    }
    int mycount_halo = 0;
    for(int k=0;k<count_halo;k+=NProcs){
        int Element = MIN(count_halo-k,NProcs);
        if(MyID < Element)
            mycount_halo ++;
    }
    int mycount = mycount_disk + mycount_bulge + mycount_halo; 
    Pall.Ntotal = Pall.NDM = mycount;
    fprintf(stderr,"mycount_disk = %d, mycount_bulge = %d, mycount_halo = %d\n",mycount_disk,mycount_bulge,mycount_halo);
    fprintf(stderr,"Ntotal = %d, NDM = %d\n",Pall.Ntotal,Pall.NDM);

    // hydro
    //// int mycount_disk_gas = mycount_disk*GasFlactionDisk; 
    //// int mycount_halo_gas = mycount_halo*GasFlactionHalo; 
    //// Pall.Nhydro = mycount_disk_gas+mycount_halo_gas;
    //// int LocalNhydro,GlobalNhydro;
    //// LocalNhydro = Pall.Nhydro;
    //// MPI_Allreduce(&LocalNhydro,&GlobalNhydro,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    //// Pall.Nhydro_t = GlobalNhydro;

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);

    GenerateStructPbody(AllocationSize);


    //double EPS_CDM = 200*PC_CGS/Pall.UnitLength;
    double EPS_CDM = 20*PC_CGS/Pall.UnitLength;
    //double EPS_CDM = 10*PC_CGS/Pall.UnitLength;
    if(MyID == MPI_ROOT_RANK)
        fprintf(stderr,"The softening length = %g pc\n",EPS_CDM*Pall.UnitLength/PC_CGS);


    // O.K.
    for(int i=0;i<NProcs;i++){
        if(MyID == i){

            struct StructReadGalactICSBodyFlag Body;
            FILE *fp;
            FileOpen(fp,fname,"rb");
            fseek(fp,sizeof(int)*3,SEEK_SET);
            int current = 0;
            for(int k=0;k<count_disk;k++){
                fread(&Body,sizeof(struct StructReadGalactICSBodyFlag),1,fp);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = Body.Pos[0];
                    Pbody[current]->Pos[1] = Body.Pos[1];
                    Pbody[current]->Pos[2] = Body.Pos[2];

                    Pbody[current]->Vel[0] = Body.Vel[0];
                    Pbody[current]->Vel[1] = Body.Vel[1];
                    Pbody[current]->Vel[2] = Body.Vel[2];

                    Pbody[current]->Mass = Body.Mass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    if(Body.GasFlag == false){
                        Pbody[current]->Type = TypeDM;
                    } else {
                        Pbody[current]->Type = TypeHydro;
                    }                        
                    Pbody[current]->GlobalID = k;
                    current ++;
                }
            }

            for(int k=0;k<count_bulge;k++){
                fread(&Body,sizeof(struct StructReadGalactICSBodyFlag),1,fp);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = Body.Pos[0];
                    Pbody[current]->Pos[1] = Body.Pos[1];
                    Pbody[current]->Pos[2] = Body.Pos[2];

                    Pbody[current]->Vel[0] = Body.Vel[0];
                    Pbody[current]->Vel[1] = Body.Vel[1];
                    Pbody[current]->Vel[2] = Body.Vel[2];

                    Pbody[current]->Mass = Body.Mass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    Pbody[current]->Type = TypeDM;

                    Pbody[current]->GlobalID = count_disk + k;
                    current ++;
                }

            }

            for(int k=0;k<count_halo;k++){
                fread(&Body,sizeof(struct StructReadGalactICSBodyFlag),1,fp);
                if(k%NProcs == MyID){
                    Pbody[current]->Pos[0] = Body.Pos[0];
                    Pbody[current]->Pos[1] = Body.Pos[1];
                    Pbody[current]->Pos[2] = Body.Pos[2];

                    Pbody[current]->Vel[0] = Body.Vel[0];
                    Pbody[current]->Vel[1] = Body.Vel[1];
                    Pbody[current]->Vel[2] = Body.Vel[2];

                    Pbody[current]->Mass = Body.Mass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = EPS_CDM;
                    if(Body.GasFlag == false){
                        Pbody[current]->Type = TypeDM;
                    } else {
                        Pbody[current]->Type = TypeHydro;
                    }

                    Pbody[current]->GlobalID = count_disk + count_bulge + k;
                    current ++;
                }
            }
            fclose(fp);

            //dprintlmpi(current);
            //dprintlmpi(Pall.Ntotal);
            //dprintlmpi(Pall.Ntotal_t);
            assert(current == Pall.Ntotal);
            Pall.Ntotal = Pall.NDM = current;

        }
        //MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MyID == MPI_ROOT_RANK)
        fprintf(stderr,"Data Load finished\n");
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    // HydroPart. 
    int count_hydro_disk = 0; 
    int count_hydro_halo = 0; 
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeHydro){
            if(Pbody[i]->GlobalID < count_disk){

                count_hydro_disk ++;
            } else {

                count_hydro_halo ++;
            }
        }
    }
    Pall.Nhydro = count_hydro_disk + count_hydro_halo;
    long int Nhydro = Pall.Nhydro;
    long int GlobalNhydro;
    MPI_Allreduce(&Nhydro,&GlobalNhydro,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    Pall.Nhydro_t = GlobalNhydro;
    fprintf(stderr,"[%02d] hydro_disk = %d, hydro_halo = %d, Nhydro = %d, Nhydro_t = %d\n",
            MPIGetMyID(),count_hydro_disk,count_hydro_halo,Pall.Nhydro,Pall.Nhydro_t);

    GenerateStructPhydro(Pall.Nhydro);
    for(int k=0;k<Pall.Nhydro;k++){
        Phydro[k] = PhydroElements+k;
    }

    // Connect Body and Hydro locally.
    int mycount_hydro = 0;
    int mycount_disk_hydro = 0;
    int mycount_halo_hydro = 0;
    for(int k=0;k<Pall.Ntotal;k++){
        if(Pbody[k]->Type == TypeHydro){
            if(Pbody[k]->GlobalID < count_disk){
                Pbody[k]->Baryon = (void *)(Phydro[mycount_hydro]);
                Phydro[mycount_hydro]->Body = Pbody[k];
                mycount_hydro ++;
                mycount_disk_hydro ++;
            } else {
                Pbody[k]->Baryon = (void *)(Phydro[mycount_hydro]);
                Phydro[mycount_hydro]->Body = Pbody[k];
                mycount_hydro ++;
                mycount_halo_hydro ++;
            }
        }
    }
    //// for(int k=0;k<Pall.Nhydro;k++){
        //// Pbody[k] = PbodyElements+k;
        //// Phydro[k] = PhydroElements+k;
        //// Pbody[k]->Baryon = (void *)(Phydro[k]);
        //// Phydro[k]->Body = Pbody[k];
    //// }



    // mycount_disk_gas -> Gas; Change type and add extra information.
    // mycount_halo_gas -> Gas; Cahnge type and add extra information.
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->Use = ON;
        Phydro[i]->Active = ON;

        Phydro[i]->Kernel = 10*PhydroBody(i)->Eps;
        if(PhydroBody(i)->GlobalID < count_disk) {
            Phydro[i]->U = Pall.ConvertTtoU*1.e+4;
            //Phydro[i]->U = Pall.ConvertTtoU*1.e+2;
        } else {
            Phydro[i]->U = Pall.ConvertTtoU*1.e+5;
            //Phydro[i]->U = Pall.ConvertTtoU*1.e+2;
        }
        // Phydro[i]->Z = Phydro[i]->ZII = Phydro[i]->ZIa = 0.01; 
        // We here assume the half solar metallicity for the ISM at initially.
#ifdef USE_CELIB //{
        //CELibSetPrimordialMetallicity(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements);
        //CELibSetSolarMetallicity(Pbody[mycount]->Mass,PbodyHydro(mycount)->Elements);

        CELibSetMetallicityWithSolarAbundancePattern(PhydroBody(i)->Mass,Phydro[i]->Elements,
                0.5*CELibGetMetalFractionForSolarChemicalComposision());
        double MassLightElements = Phydro[i]->Elements[CELibYield_H]+Phydro[i]->Elements[CELibYield_He];
        double Z = (PhydroBody(i)->Mass-MassLightElements)/PhydroBody(i)->Mass;
        Phydro[i]->Z = Phydro[i]->ZII = Phydro[i]->ZIa = Z;
#else 
        Phydro[i]->Z   = 0.01;
        Phydro[i]->ZII = 0.01;
        Phydro[i]->ZIa = 0.00;
#endif // USE_CELIB //}

#if (UseSFModelSpawn) 
        Phydro[i]->SpawnMass = PhydroBody(i)->Mass/(double)MaxSpawnTimes;
#endif

        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
    }

#if 0
    for(int i=0;i<mycount_disk_hydro+mycount_halo_hydro;i++){
        //PhydroBody(i)->Type = TypeHydro;

        Phydro[i]->Use = ON;
        Phydro[i]->Active = ON;
        if(i<mycount_disk_hydro){
            Phydro[i]->U = Pall.ConvertTtoU*1.e+4;
        } else {
            Phydro[i]->U = Pall.ConvertTtoU*1.e+5;
        }
        Phydro[i]->Kernel = 10*PhydroBody(i)->Eps;
#if (UseSFModelSpawn) 
        Phydro[i]->SpawnMass = PhydroBody(i)->Mass/(double)MaxSpawnTimes;
#endif
        Phydro[i]->Z = Phydro[i]->ZII = Phydro[i]->ZIa = 0.01; 
        // We here assume the half solar metallicity for the ISM at initially.

        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
    }
#endif

    fprintf(stderr,"Thalo = %g K.\n",
            (0.59*PROTON_MASS_CGS/(3.0*BOLTZMANN_CONSTANT_CGS))*(GRAVITY_CONSTANT_CGS*1.2e+11*MSUN_CGS/(100*KPC_CGS)));

    ReConnectPointers();
    UpdateTotalNumber();


    //MPI_Barrier(MPI_COMM_WORLD);
    //exit(1);

    Pall.RunStatus = NewSimulation;
    // hydro parameters
    // Pall.Ns = 32;
    // Pall.Npm = 2;
    Pall.Ns = 128;
    Pall.Npm = 8;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    // Pall.hubble = 0.5;
    // Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    // Pall.OmegaM = 1.0;
    // Pall.OmegaL = 0.0;
    // Pall.OmegaB = 0.0;
    Pall.InitialRedshift = Pall.Redshift = 0;

    Pall.TEnd = 1.5*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g\n",Pall.TCurrent,Pall.TEnd);

    Pall.AdaptiveSofteningFactor = 1.e0;
    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = 0.1*GIGAYEAR_CGS/Pall.UnitTime; 
    //Pall.OutPutInterval = 0.01*GIGAYEAR_CGS/Pall.UnitTime; 
    Pall.OutPutInterval = Pall.TEnd/1500.0; 
    //Pall.OutPutInterval = Pall.TEnd/150.0; 

    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/ID.ASCII");
    strcpy(Pall.BaseFileName,"./data/ID");
    strcpy(Pall.RestartFileName,"./data/ID.dump");


    return;
}

struct StructDM{
    unsigned long int   GlobalID;
    float Mass;
    float Pos[3];
    float Vel[3];
    float Eps;
};

struct StructGas{
    unsigned long int   GlobalID;
    float Mass;
    float Pos[3];
    float Vel[3];
    float Eps;
    float Kernel;
    float U;
    float Z;
    float ZII;
    float ZIa;
    short   SpawnTimes;
    double  SpawnMass;
};

struct StructStars{
    unsigned long int   GlobalID;
    float Mass;
    float Pos[3];
    float Vel[3];
    float Eps;


    short   IMFTYPE;
    bool    TypeII;
    bool    TypeIa;
    short   NthChildren;
    unsigned long int   ParentGlobalID;
    double  InitialMass;    // Initial Mass 
    double  sMass;           // Current Mass
    double  FormationTime;
    double  Z;              // Metallicity of gas.
    double  ZII;            // The mass of metal by TypeII.
    double  ZIa;            // The mass of metal by TypeIa.
    double  TempForm;       // Temperature(t=FormationTime)
    double  RhoForm;        // Density(t=FormationTime) // UnitMass/UnitLength^3
};



void ReadMergerRestart(char fname[]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = KPC_CGS;
    Pall.UnitMass = MSUN_CGS;
    Pall.UnitTime = KPC_CGS/(VELOCITY_KMS_CGS);
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    // Set cosmological parameters!
    
    // read body! 
    int count_hydro,count_stars,count_dm;
    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,fname,"rb");
        fread(&count_hydro,sizeof(int),1,fp);
        fread(&count_stars,sizeof(int),1,fp);
        fread(&count_dm,sizeof(int),1,fp);
        fclose(fp);
    }
    MPI_Bcast(&count_hydro,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&count_stars,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&count_dm,1,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);

    fprintf(stderr,"[%02d] Nhydro,Nstars,Ndm = %d %d %d\n",MyID,count_hydro,count_stars,count_dm);

    Pall.Nhydro_t = count_hydro;
    Pall.Nstars_t = count_stars;
    Pall.NDM_t = count_dm;
    Pall.Ntotal_t = count_hydro + count_stars + count_dm;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"Ntotal_t = %ld, NDM_t = %ld\n",Pall.Ntotal_t,Pall.NDM_t);

    int mycount_hydro = 0;
    for(int k=0;k<count_hydro;k+=NProcs){
        int Element = MIN(count_hydro-k,NProcs);
        if(MyID < Element)
            mycount_hydro ++;
    }
    int mycount_stars = 0;
    for(int k=0;k<count_stars;k+=NProcs){
        int Element = MIN(count_stars-k,NProcs);
        if(MyID < Element)
            mycount_stars ++;
    }
    int mycount_dm = 0;
    for(int k=0;k<count_dm;k+=NProcs){
        int Element = MIN(count_dm-k,NProcs);
        if(MyID < Element)
            mycount_dm ++;
    }

    int mycount = mycount_hydro + mycount_stars + mycount_dm; 
    Pall.Nhydro = mycount_hydro;
    Pall.Nstars = mycount_stars;
    Pall.NDM = mycount_dm;
    Pall.Ntotal = mycount;
    fprintf(stderr,"mycount_hydro = %d, mycount_stars = %d, mycount_dm = %d\n",mycount_hydro,mycount_stars,mycount_dm);
    fprintf(stderr,"Ntotal = %ld, NDM = %ld\n",Pall.Ntotal,Pall.NDM);

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro((int)(ForAngelsShare*mycount_hydro));
    GenerateStructPstar((int)(ForAngelsShare*mycount_stars));

    for(int i=0;i<mycount_hydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    for(int i=0;i<mycount_stars;i++){
        Pbody[i+mycount_hydro]->Baryon = (void *)(Pstar[i]);
        Pstar[i]->Body = Pbody[i+mycount_hydro];
    }

    for(int i=0;i<NProcs;i++){
        if(MyID == i){
            FILE *fp;
            FileOpen(fp,fname,"rb");
            fseek(fp,sizeof(int)*3,SEEK_SET);
            int current = 0;

            struct StructGas Gas;
            for(int k=0;k<count_hydro;k++){
                fread(&Gas,sizeof(struct StructGas),1,fp);
                if(k%NProcs == MyID){
                    Pbody[current]->GlobalID = Gas.GlobalID;
                    Pbody[current]->Pos[0] = Gas.Pos[0];
                    Pbody[current]->Pos[1] = Gas.Pos[1];
                    Pbody[current]->Pos[2] = Gas.Pos[2];

                    Pbody[current]->Vel[0] = Gas.Vel[0];
                    Pbody[current]->Vel[1] = Gas.Vel[1];
                    Pbody[current]->Vel[2] = Gas.Vel[2];

                    Pbody[current]->Mass = Gas.Mass;
                    PbodyHydro(current)->Kernel = Gas.Kernel;
                    PbodyHydro(current)->U = Gas.U;
                    PbodyHydro(current)->Z = Gas.Z;
                    PbodyHydro(current)->ZII = Gas.ZII;
                    PbodyHydro(current)->ZIa = Gas.ZIa;
                    PbodyHydro(current)->SpawnTimes = Gas.SpawnTimes;
                    PbodyHydro(current)->SpawnMass = Gas.SpawnMass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = Gas.Eps;
                    Pbody[current]->Type = TypeHydro;
                    PbodyHydro(current)->Use = ON;
                    PbodyHydro(current)->Active = ON;
                    current ++;
                }
            }

            struct StructStars Stars;
            for(int k=0;k<count_stars;k++){
                fread(&Stars,sizeof(struct StructStars),1,fp);
                if(k%NProcs == MyID){
                    Pbody[current]->GlobalID = Stars.GlobalID;
                    Pbody[current]->Pos[0] = Stars.Pos[0];
                    Pbody[current]->Pos[1] = Stars.Pos[1];
                    Pbody[current]->Pos[2] = Stars.Pos[2];

                    Pbody[current]->Vel[0] = Stars.Vel[0];
                    Pbody[current]->Vel[1] = Stars.Vel[1];
                    Pbody[current]->Vel[2] = Stars.Vel[2];
                    Pbody[current]->Mass = Stars.Mass;

                    PbodyStar(current)->IMFTYPE = Stars.IMFTYPE;
                    PbodyStar(current)->TypeII = Stars.TypeII;
                    PbodyStar(current)->TypeIa = Stars.TypeIa;
                    PbodyStar(current)->NthChildren = Stars.NthChildren;
                    PbodyStar(current)->ParentGlobalID = Stars.ParentGlobalID;
                    PbodyStar(current)->InitialMass = Stars.InitialMass;
                    PbodyStar(current)->Mass = Stars.sMass;
                    PbodyStar(current)->FormationTime = Stars.FormationTime;
                    PbodyStar(current)->Z = Stars.Z;
                    PbodyStar(current)->ZII = Stars.ZII;
                    PbodyStar(current)->ZIa = Stars.ZIa;
                    PbodyStar(current)->TempForm = Stars.TempForm;
                    PbodyStar(current)->RhoForm = Stars.RhoForm;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = Stars.Eps;
                    Pbody[current]->Type = TypeStar;

                    PbodyStar(current)->Use = ON;

                    current ++;
                }
            }

            struct StructDM DM;
            for(int k=0;k<count_dm;k++){
                fread(&DM,sizeof(struct StructDM),1,fp);
                if(k%NProcs == MyID){
                    Pbody[current]->GlobalID = DM.GlobalID;
                    Pbody[current]->Pos[0] = DM.Pos[0];
                    Pbody[current]->Pos[1] = DM.Pos[1];
                    Pbody[current]->Pos[2] = DM.Pos[2];

                    Pbody[current]->Vel[0] = DM.Vel[0];
                    Pbody[current]->Vel[1] = DM.Vel[1];
                    Pbody[current]->Vel[2] = DM.Vel[2];

                    Pbody[current]->Mass = DM.Mass;

                    Pbody[current]->Active = ON;
                    Pbody[current]->Use = ON;
                    Pbody[current]->Eps = DM.Eps;
                    Pbody[current]->Type = TypeDM;
                    current ++;
                }
            }
            fclose(fp);

            assert(current == Pall.Ntotal);
            Pall.Ntotal = current;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if(MyID == MPI_ROOT_RANK)
        fprintf(stderr,"Data Load finished\n");
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    //fprintf(stderr,"%ld %ld %ld\n",Pbody[0]->GlobalID,PhydroBody(0)->GlobalID,PbodyElements[0].GlobalID);
    //fprintf(stderr,"%d %d %d\n",Pbody[0]->Use,Phydro[0]->Use,PhydroElements[0].Use);

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
    }

    ReConnectPointers();
    UpdateTotalNumber();

    Pall.RunStatus = NewSimulation;
    // hydro parameters
    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    // Pall.hubble = 0.5;
    // Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    // Pall.OmegaM = 1.0;
    // Pall.OmegaL = 0.0;
    // Pall.OmegaB = 0.0;
    Pall.InitialRedshift = Pall.Redshift = 0;

    Pall.TEnd = 0.01*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g\n",Pall.TCurrent,Pall.TEnd);

    Pall.AdaptiveSofteningFactor = 1.e0;
    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = 0.1*GIGAYEAR_CGS/Pall.UnitTime; 
    //Pall.OutPutInterval = 0.01*GIGAYEAR_CGS/Pall.UnitTime; 
    //Pall.OutPutInterval = Pall.TEnd/1500.0; 
    Pall.OutPutInterval = Pall.TEnd/1.0; 

    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/ID.ASCII");
    strcpy(Pall.BaseFileName,"./data/ID");
    strcpy(Pall.RestartFileName,"./data/ID.dump");


    return;
}
#endif

#ifdef TASK_DICE_RUN //{

// #define FOR_DWARF

struct StructDiceRun{
    int SimulationModel; // 0 = MW, 1 = Dwarf
    int MetallicityType; // 0 = Solar(MW), 1 = 0.1 solar(Dwarf), 2 = Primordial, 3 = follow the init condition
    int UseEstimateXFe;  // 0 = not use, 1 = use, sofar only available for the case with MetallicityType=3
    double Eps[4]; // Eps for particles of gas, dm, (halo) stars, bulge stars.
} DiceRunParameters;

static void LoadDiceRunParameters(void){

    char fname[MaxCharactersInLine] = "Run.params";

    if(CheckFile(fname)){ // load
        FILE *fp;
        FileOpen(fp,fname,"r");
        fscanf(fp,"%d",&DiceRunParameters.SimulationModel);
        fscanf(fp,"%d",&DiceRunParameters.MetallicityType);
        fscanf(fp,"%le",DiceRunParameters.Eps);
        fscanf(fp,"%le",DiceRunParameters.Eps+1);
        fscanf(fp,"%le",DiceRunParameters.Eps+2);
        fscanf(fp,"%le",DiceRunParameters.Eps+3);
        fclose(fp);

        for(int i=0;i<4;i++){
            DiceRunParameters.Eps[i] *= PC_CGS/Pall.UnitLength;
        }

        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"\033[1;33m== Reload parameters (use Run.params) ! ==\033[0m\n");
            fprintf(stderr,"\033[1;33m== SimulationModel %d\033[0m\n",DiceRunParameters.SimulationModel);
            fprintf(stderr,"\033[1;33m== MetallicityType %d\033[0m\n",DiceRunParameters.MetallicityType);
            fprintf(stderr,"\033[1;33m== Eps[0] (for gas) %g\033[0m\n",DiceRunParameters.Eps[0]);
            fprintf(stderr,"\033[1;33m== Eps[1] (for dm) %g\033[0m\n",DiceRunParameters.Eps[1]);
            fprintf(stderr,"\033[1;33m== Eps[2] (for (halo) stars) %g\033[0m\n",DiceRunParameters.Eps[2]);
            fprintf(stderr,"\033[1;33m== Eps[3] (for bulge stars) %g\033[0m\n",DiceRunParameters.Eps[3]);
        }
    } else { // set defaults
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"\033[1;33m== Use pre-existing parameters (not use Run.params) ! ==\033[0m\n");
            fprintf(stderr,"\033[1;33m== SimulationModel %d\033[0m\n",DiceRunParameters.SimulationModel);
            fprintf(stderr,"\033[1;33m== MetallicityType %d\033[0m\n",DiceRunParameters.MetallicityType);
            fprintf(stderr,"\033[1;33m== Eps[0] (for gas) %g\033[0m\n",DiceRunParameters.Eps[0]);
            fprintf(stderr,"\033[1;33m== Eps[1] (for dm) %g\033[0m\n",DiceRunParameters.Eps[1]);
            fprintf(stderr,"\033[1;33m== Eps[2] (for (halo) stars) %g\033[0m\n",DiceRunParameters.Eps[2]);
            fprintf(stderr,"\033[1;33m== Eps[3] (for bulge stars) %g\033[0m\n",DiceRunParameters.Eps[3]);
        }
    }

    return ;
}

#define USE_DICE_EPS_FOLLOW_MASS
#define USE_DICE_BARYON_SOFTENING

static double ComputeEps(const double Eps, const double Mass, const int Type){

#ifdef USE_DICE_EPS_FOLLOW_MASS //{

#define Mvir (1.e+12)
#ifdef USE_DICE_BARYON_SOFTENING //{
    if(Type == TypeHydro){
        // 100*mass/eps^3 = rhoth;
        // eps = (100*mass/rhoth)^{1/3}
        // need to consider 4pi/3 ?
        double rhoth = SFCONDITION_DENSITY_CRITERION/Pall.ConvertNumberDensityToCGS;
        return cbrt(100*Mass/rhoth);
    } else {
        double Mass_in_Solar = Mass*Pall.UnitMass/MSUN_CGS;
        double Eps_new = 30*sqrt(Mass_in_Solar/1.e+3)*pow((Mvir/1.e+12),0.2)*PC_CGS/Pall.UnitLength;
        return Eps_new/3.0;
    }
#else  // USE_DICE_BARYON_SOFTENING //}//{
    double Mass_in_Solar = Mass*Pall.UnitMass/MSUN_CGS;
    double Eps_new = 30*sqrt(Mass_in_Solar/1.e+3)*pow((Mvir/1.e+12),0.2)*PC_CGS/Pall.UnitLength;
    return Eps_new/3.0;
#endif // USE_DICE_BARYON_SOFTENING //}


#else // USE_DICE_EPS_FOLLOW_MASS //}//{
    return Eps;
#endif  // USE_DICE_EPS_FOLLOW_MASS //}
#undef Mvir

}

static void InsertDiceData(const int TargetIndex, float Pos[], float Vel[], float Mass, float Rho, float U,
        const double Age, const double Metal, const int Type, const int GlobalID, const double Eps){

    Pbody[TargetIndex]->GlobalID = GlobalID;
    Pbody[TargetIndex]->Pos[0] = Pos[0];
    Pbody[TargetIndex]->Pos[1] = Pos[1];
    Pbody[TargetIndex]->Pos[2] = Pos[2];

    Pbody[TargetIndex]->Vel[0] = Vel[0];
    Pbody[TargetIndex]->Vel[1] = Vel[1];
    Pbody[TargetIndex]->Vel[2] = Vel[2];
    
    Pbody[TargetIndex]->Mass = Mass;

    /*
    fprintf(stderr,"%g | %g %g %g | %g %g %g \n",Mass,Pos[0],Pos[1],Pos[2],
            Vel[0],Vel[1],Vel[2]);
    */

    Pbody[TargetIndex]->Active = ON;
    Pbody[TargetIndex]->Use = ON;
    Pbody[TargetIndex]->Type = Type;
    //Pbody[TargetIndex]->Eps = Eps;
    Pbody[TargetIndex]->Eps = ComputeEps(Eps,Mass,Type);

    if(Type == TypeHydro){
        PbodyHydro(TargetIndex)->Use = ON;
        PbodyHydro(TargetIndex)->Rho = Rho;
        PbodyHydro(TargetIndex)->U = U;
        PbodyHydro(TargetIndex)->Kernel = 2.0*Eps;
        PbodyHydro(TargetIndex)->Alpha = 1.e0;
#if (UseSFModelSpawn) 
        PbodyHydro(TargetIndex)->SpawnMass = Pbody[TargetIndex]->Mass/(double)MaxSpawnTimes;
#endif
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
        PbodyHydro(TargetIndex)->G0 = 1;
        PbodyHydro(TargetIndex)->fH2 = 0.1;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

#ifdef USE_CELIB //{
#if 1
        if(DiceRunParameters.MetallicityType == 0){ // for mw
            CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements,
                 0.1*CELibGetMetalFractionForSolarChemicalComposision());
        }else if(DiceRunParameters.MetallicityType == 1){ // for dwarf
            CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements,
                 0.1*CELibGetMetalFractionForSolarChemicalComposision());
        }else if(DiceRunParameters.MetallicityType == 2){ // for primordial runs
            CELibSetPrimordialMetallicity(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements);
        }else if(DiceRunParameters.MetallicityType == 3){ // follow the init file
            if(DiceRunParameters.UseEstimateXFe == 0){
                CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements,Metal);
            } else if(DiceRunParameters.UseEstimateXFe == 1){
                double FeH = Metal/0.0134;
                EstimateAbudancePatternForMass(PbodyHydro(TargetIndex)->Elements, FeH,Pbody[TargetIndex]->Mass,XFeFeH_SAGACompact);
#if 0
                {
                    double sum = 0.e0;
                    for(int l=0;l<13;l++){
                        sum += PbodyHydro(TargetIndex)->Elements[l];
                    }
                    // fprintf(stderr,"-- %g %g %g %g\n",Pbody[TargetIndex]->Mass,sum,PbodyHydro(TargetIndex)->Elements[2],FeH);
                }
#endif
            }
        }
#else //
        //CELibSetPrimordialMetallicity(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements);
#ifdef FOR_DWARF
        CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements,
                 0.1*CELibGetMetalFractionForSolarChemicalComposision());
#else
        //CELibSetSolarMetallicity(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements);
        //CELibSetPrimordialMetallicity(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements);
        //CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements,
                 //0.01*CELibGetMetalFractionForSolarChemicalComposision());
#endif
        CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyHydro(TargetIndex)->Elements,Metal);
#endif




        double MassLightElements = PbodyHydro(TargetIndex)->Elements[CELibYield_H]
                                  +PbodyHydro(TargetIndex)->Elements[CELibYield_He];
        double Z = (Pbody[TargetIndex]->Mass-MassLightElements)/Pbody[TargetIndex]->Mass;
        PbodyHydro(TargetIndex)->Z = PbodyHydro(TargetIndex)->ZII = PbodyHydro(TargetIndex)->ZIa = Z;
#else 
        PbodyHydro(TargetIndex)->Z   = 0.02;
        PbodyHydro(TargetIndex)->ZII = 0.02;
        PbodyHydro(TargetIndex)->ZIa = 0.00;
#endif // USE_CELIB //}
    } else if(Type == TypeStar){
#include "StellarFeedback.h"
        PbodyStar(TargetIndex)->Use = ON;
        PbodyStar(TargetIndex)->IMFTYPE = IMFTYPE_SP;
        PbodyStar(TargetIndex)->Mass = Pbody[TargetIndex]->Mass;
        PbodyStar(TargetIndex)->InitialMass = Pbody[TargetIndex]->Mass;
        PbodyStar(TargetIndex)->FormationTime = Age;
        //CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyStar(TargetIndex)->Elements,Metal);
        if(DiceRunParameters.UseEstimateXFe == 0){
            CELibSetMetallicityWithSolarAbundancePattern(Pbody[TargetIndex]->Mass,PbodyStar(TargetIndex)->Elements,Metal);
        } else if(DiceRunParameters.UseEstimateXFe == 1){
            double FeH = Metal/0.0134;
            EstimateAbudancePatternForMass(PbodyStar(TargetIndex)->Elements,FeH,Pbody[TargetIndex]->Mass,XFeFeH_SAGACompact);
#if 0
            {
                double sum = 0.e0;
                for(int l=0;l<13;l++){
                    sum += PbodyStar(TargetIndex)->Elements[l];
                }
                double Znew = (sum-PbodyStar(TargetIndex)->Elements[0]-PbodyStar(TargetIndex)->Elements[1])/sum;
                //fprintf(stderr,"++ %g %g %g %g\n",Pbody[TargetIndex]->Mass,sum,PbodyStar(TargetIndex)->Elements[2],FeH);
                double ElementsDummy[CELibYield_Number] = {0.e0};
                CELibSetSolarMetallicity(Pbody[TargetIndex]->Mass,ElementsDummy);
                fprintf(stderr,"%g %g %g %g %g %g %g %g %g %g %g\n",
                        NORM(Pbody[TargetIndex]->Pos),Pbody[TargetIndex]->Mass,sum,PbodyStar(TargetIndex)->Elements[2],FeH,Znew,
                        PbodyStar(TargetIndex)->Elements[CELibYield_Fe],
                        ElementsDummy[CELibYield_Fe],PbodyStar(TargetIndex)->Elements[CELibYield_Fe]/ElementsDummy[CELibYield_Fe],
                        PbodyStar(TargetIndex)->Elements[CELibYield_Fe]/sum,
                        ElementsDummy[CELibYield_Fe]/sum);
                /*
                        */
            }
#endif
        }

        PbodyStar(TargetIndex)->Z = (PbodyStar(TargetIndex)->Mass-PbodyStar(TargetIndex)->Elements[CELibYield_H]-PbodyStar(TargetIndex)->Elements[CELibYield_He])/PbodyStar(TargetIndex)->Mass;
        PbodyStar(TargetIndex)->ZII = PbodyStar(TargetIndex)->Z;
        PbodyStar(TargetIndex)->ZIa = 0;

        if(fabs(Age)*Pall.UnitTime < 40.0*MEGAYEAR_CGS){
            PbodyStar(TargetIndex)->TypeII = false;
        } else {
            PbodyStar(TargetIndex)->TypeII = true;
        }
        PbodyStar(TargetIndex)->TypeIa = false;
        PbodyStar(TargetIndex)->NthChildren = 0;
        PbodyStar(TargetIndex)->ParentGlobalID = GlobalID;

#ifdef USE_CELIB
        PbodyStar(TargetIndex)->IMFTYPE = StellarFeedbackGetIMFType();
        PbodyStar(TargetIndex)->SNIaCount = -1;
#if 1
        if(PbodyStar(TargetIndex)->TypeII == true){
#if 0
            do{
                PbodyStar(TargetIndex)->SNIaCount ++;
                PbodyStar(TargetIndex)->EventTime = PbodyStar(TargetIndex)->FormationTime 
                                   +StellarFeedbackGetNextEventTime(CELib_Feedback_TypeIa,PbodyStar(TargetIndex)->Z,PbodyStar(TargetIndex)->InitialMass,0)
                                    *YEAR_CGS/Pall.UnitTime;
            } while(PbodyStar(TargetIndex)->EventTime < 0);
#endif
            PbodyStar(TargetIndex)->SNIICount = NONE;
            PbodyStar(TargetIndex)->EventTimeSNII = 100*GIGAYEAR_CGS/Pall.UnitTime;
            PbodyStar(TargetIndex)->SNIaCount = 100000;
            PbodyStar(TargetIndex)->EventTimeSNIa = 100*GIGAYEAR_CGS/Pall.UnitTime;
        } else {
            PbodyStar(TargetIndex)->SNIICount = 0;
            PbodyStar(TargetIndex)->SNIaCount = 0;
            PbodyStar(TargetIndex)->EventTimeSNII = PbodyStar(TargetIndex)->FormationTime 
                            +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                                .R = gsl_rng_uniform(RandomGenerator),
                                .InitialMass_in_Msun = PbodyStar(TargetIndex)->InitialMass*Pall.UnitMass/MSUN_CGS,
                                .Metallicity = PbodyStar(TargetIndex)->Z,
                                .Count = 0,
                                },CELibFeedbackType_SNII)
                                *YEAR_CGS/Pall.UnitTime;
                               //+StellarFeedbackGetNextEventTime(CELibFeedbackType_SNII,PbodyStar(TargetIndex)->Z,PbodyStar(TargetIndex)->InitialMass,0)
                                //*YEAR_CGS/Pall.UnitTime;
        }
#ifdef USE_CELIB_AGB //{
#if 0
        do{
            PbodyStar(TargetIndex)->AGBCount ++;
            PbodyStar(TargetIndex)->EventTimeAGB = PbodyStar(TargetIndex)->FormationTime 
                               +StellarFeedbackGetNextEventTime(CELib_Feedback_AGB,PbodyStar(TargetIndex)->Z,PbodyStar(TargetIndex)->InitialMass,0)
                                *YEAR_CGS/Pall.UnitTime;
        } while(PbodyStar(TargetIndex)->EventTimeAGB < 0);
#endif
        PbodyStar(TargetIndex)->EventTimeAGB = 100*GIGAYEAR_CGS/Pall.UnitTime;
        PbodyStar(TargetIndex)->AGBCount = 100000000;
#endif //USE_CELIB_AGB //}
#else
    if(PbodyStar(TargetIndex)->TypeII == true){
        PbodyStar(TargetIndex)->SNIaCount = 1.e9;
        PbodyStar(TargetIndex)->EventTime = 100*GIGAYEAR_CGS/Pall.UnitTime;
    } else {
        PbodyStar(TargetIndex)->SNIaCount = -1;
        PbodyStar(TargetIndex)->EventTime = PbodyStar(TargetIndex)->FormationTime 
                            +CELibGetNextEventTime((struct CELibStructNextEventTimeInput){
                                .R = gsl_rng_uniform(RandomGenerator),
                                .InitialMass_in_Msun = PbodyStar(TargetIndex)->InitialMass*Pall.UnitMass/MSUN_CGS,
                                .Metallicity = PbodyStar(TargetIndex)->Z,
                                .Count = 0,
                                },CELibFeedbackType_SNII)
                                *YEAR_CGS/Pall.UnitTime;
                           //+StellarFeedbackGetNextEventTime(CELibFeedbackType_SNII,PbodyStar(TargetIndex)->Z,PbodyStar(TargetIndex)->InitialMass,0)
                            //*YEAR_CGS/Pall.UnitTime;
    }
    PbodyStar(TargetIndex)->EventTimeAGB = 100*GIGAYEAR_CGS/Pall.UnitTime;
    PbodyStar(TargetIndex)->AGBCount = 1.e9;
#endif
#endif //USE_CELIB_AGB //}

#ifdef USE_CELIB_NSM //{
        PbodyStar(TargetIndex)->EventTimeNSM = 100*GIGAYEAR_CGS/Pall.UnitTime;
        PbodyStar(TargetIndex)->NSMCount = 100000000;
#endif // USE_CELIB_NSM //}

    }

    return ;
}

static void ReadDiceLine(char fname[], const int Nlines, const int mycount_comp[], const double Eps[]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int Number = mycount_comp[0] + mycount_comp[1] + mycount_comp[2] + mycount_comp[3];
    dprintlmpi(Nlines);
    dprintlmpi(Number);

    FILE *fp;
    FileOpen(fp,fname,"r");
    fscanf(fp,"%*d %*d %*d %*d");
    fscanf(fp,"%*e %*e %*e %*e");

    int counter_total = 0;
    int counter_hydro = 0;
    int counter_dm = 0;
    int counter_stars = 0;
    for(int i=0;i<Nlines;i++){
        float Pos[3],Vel[3],Mass,U,Rho,Metal,Age;
        int Type;
        fscanf(fp,"%e %e %e %e %e %e %e %e %e %d %e %e",
                Pos,Pos+1,Pos+2,Vel,Vel+1,Vel+2,
                &U,&Rho,&Mass,&Type,&Metal,&Age);
#if 1
        Vel[0] = Vel[0]*1.e5/(Pall.UnitLength/Pall.UnitTime); // km/s -> cm/s -> L/T
        Vel[1] = Vel[1]*1.e5/(Pall.UnitLength/Pall.UnitTime);
        Vel[2] = Vel[2]*1.e5/(Pall.UnitLength/Pall.UnitTime);
#endif

#ifdef FOR_DWARF
        U = Pall.ConvertTtoU*3.e+2;
        Metal = 0.01304*0.1;
#endif

        if(i%NProcs == MyID){
            if(Type == 0){ // hydro
                InsertDiceData(counter_hydro,Pos,Vel,Mass,Rho,U,Age,Metal,TypeHydro,i,Eps[Type]);
                counter_hydro ++;
            } else if(Type == 1){ //DM
                InsertDiceData(counter_dm+mycount_comp[0]+mycount_comp[2]+mycount_comp[3],Pos,Vel,Mass,Rho,U,Age,Metal,TypeDM,i,Eps[Type]);
                counter_dm ++;
            } else {
                InsertDiceData(counter_stars+mycount_comp[0],Pos,Vel,Mass,Rho,U,Age,Metal,TypeStar,i,Eps[Type]);
                counter_stars ++;
            }
            counter_total ++;
        }
    }
    fclose(fp);

    dprintlmpi(counter_hydro);
    dprintlmpi(counter_dm);
    dprintlmpi(counter_stars);
    assert(Number == counter_total);
    // exit(1);

    return ;
}

static void ReadDiceLineBinary(char fname[], const int Nlines, const int mycount_comp[], const double Eps[]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    int Number = mycount_comp[0] + mycount_comp[1] + mycount_comp[2] + mycount_comp[3];
    // dprintlmpi(Nlines);
    // dprintlmpi(Number);

    int Numbers[4] = {0,0,0,0};
    double Units[4] = {0.e0,0.e0,0.e0,0.e0};

    FILE *fp;
    FileOpen(fp,fname,"r");
    fread(Numbers+0,sizeof(int),1,fp);
    fread(Numbers+1,sizeof(int),1,fp);
    fread(Numbers+2,sizeof(int),1,fp);
    fread(Numbers+3,sizeof(int),1,fp);
    fread(Units+0,sizeof(double),1,fp);
    fread(Units+1,sizeof(double),1,fp);
    fread(Units+2,sizeof(double),1,fp);
    fread(Units+3,sizeof(double),1,fp);

    int counter_total = 0;
    int counter_hydro = 0;
    int counter_dm = 0;
    int counter_stars = 0;
    for(int i=0;i<Nlines;i++){
        float Pos[3],Vel[3],Mass,U,Rho,Metal,Age;
        int Type;

        fread(Pos,sizeof(float),3,fp);
        fread(Vel,sizeof(float),3,fp);
        fread(&U,sizeof(float),1,fp);
        fread(&Rho,sizeof(float),1,fp);
        fread(&Mass,sizeof(float),1,fp);
        fread(&Type,sizeof(int),1,fp);
        fread(&Metal,sizeof(float),1,fp);
        fread(&Age,sizeof(float),1,fp);

#if 1
        Vel[0] = Vel[0]*1.e5/(Pall.UnitLength/Pall.UnitTime); // km/s -> cm/s -> L/T
        Vel[1] = Vel[1]*1.e5/(Pall.UnitLength/Pall.UnitTime);
        Vel[2] = Vel[2]*1.e5/(Pall.UnitLength/Pall.UnitTime);
#endif

        // U = Pall.ConvertTtoU*3.e+2;
        // Metal = 0.01304*0.1;
        if(i%NProcs == MyID){
            if(Type == 0){ // hydro
                InsertDiceData(counter_hydro,Pos,Vel,Mass,Rho,U,Age,Metal,TypeHydro,i,Eps[Type]);
                counter_hydro ++;
            } else if(Type == 1){ //DM
                InsertDiceData(counter_dm+mycount_comp[0]+mycount_comp[2]+mycount_comp[3],Pos,Vel,Mass,Rho,U,Age,Metal,TypeDM,i,Eps[Type]);
                counter_dm ++;
            } else {
                InsertDiceData(counter_stars+mycount_comp[0],Pos,Vel,Mass,Rho,U,Age,Metal,TypeStar,i,Eps[Type]);
                counter_stars ++;
            }
            counter_total ++;
        }
    }
    fclose(fp);

    dprintlmpi(counter_hydro);
    dprintlmpi(counter_dm);
    dprintlmpi(counter_stars);
    assert(Number == counter_total);
    // exit(1);

    return ;
}

void NumberCounter(int Numbers[], char fname[]){

    for(int i=0;i<8;i++)
        Numbers[i] = 0;

    FILE *fp;
    FileOpen(fp,fname,"r");
    int Lines[4];
    fscanf(fp,"%d %d %d %d",Lines,Lines+1,Lines+2,Lines+3);
    int Nlines = Lines[0]+Lines[1]+Lines[2]+Lines[3];
    fscanf(fp,"%*e %*e %*e %*e");
    int counter = 0;
    int OldType = 0;
    for(int i=0;i<Nlines;i++){
        int Type;
        fscanf(fp,"%*e %*e %*e %*e %*e %*e %*e %*e %*e %d %*e %*e",&Type);
        if(i>0){
            if(OldType != Type)
                counter ++;
            Numbers[counter] ++;
        } else {
            Numbers[counter] ++;
        }
        OldType = Type;
    }
    fclose(fp);
    for(int i=0;i<8;i++)
        dprintlmpi(Numbers[i]);

    fflush(NULL);
    exit(1);
    return ;
}

void ReadDiceASCII(char fname[]){

#ifdef USE_CELIB
    InitializeStellarFeedback();
#endif

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    
    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    int Numbers[4] = {0,0,0,0};
    double Units[4] = {0.e0,0.e0,0.e0,0.e0};
    FILE *fp;
    FileOpen(fp,fname,"r");
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fscanf(fp,"%d %d %d %d",Numbers,Numbers+1,Numbers+2,Numbers+3);
        fscanf(fp,"%le %le %le %le",Units,Units+1,Units+2,Units+3);
    }
    fclose(fp);
    MPI_Bcast(&Numbers,4,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&Units,4,MPI_DOUBLE,MPI_ROOT_RANK,MPI_COMM_WORLD);

    // Set Shuffule/Lock information
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Particle numbers, Ngas,Ndisk,Nhalo,Nspheroid = %d %d %d %d\n",
                Numbers[0],Numbers[1],Numbers[2],Numbers[3]);
        fprintf(stderr,"Simulation units, Mass,Length,Time = %g [g],  %g [cm], %g [s]\n",
                Units[0],Units[1],Units[2]);
    }
    Pall.Ntotal_t = Numbers[0]+Numbers[1]+Numbers[2]+Numbers[3];


    // Pall.UnitLength = KPC_CGS;
    // Pall.UnitMass = MSUN_CGS;
    // Pall.UnitTime = KPC_CGS/(VELOCITY_KMS_CGS);
    Pall.UnitMass = Units[0];
    Pall.UnitLength = Units[1];
    Pall.UnitTime = Units[2];
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    //SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
    SetVariableViscosityParameters(VARIABLE_VISCOSITY_MIN,VARIABLE_VISCOSITY_MAX,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

#if 0
    int mycount_comp[4] = {0,0,0,0};
    int mycount = 0;
    for(int i=0;i<4;i++){
        for(int k=0;k<Numbers[i];k+=NProcs){
            int Element = MIN(Numbers[i]-k,NProcs);
            if(MyID < Element)
                mycount_comp[i] ++;
        }
        mycount += mycount_comp[i]; 
    }
#else
    int mycount_comp[4] = {0,0,0,0};
    int mycount = 0;

    // int NumComp[8];
    // if(MPIGetMyID() == MPI_ROOT_RANK){
        // NumberCounter(NumComp,fname);
    // }

#if 0
    for(int i=0;i<4;i++){
        for(int k=0;k<Numbers[i];k+=NProcs){
            int Element = MIN(Numbers[i]-k,NProcs);
            if(MyID < Element)
                mycount_comp[i] ++;
        }
        mycount += mycount_comp[i]; 
    }
#endif
    for(int k=0;k<Pall.Ntotal_t;k++){
        if(k%NProcs == MyID){
            if(k<Numbers[0]){
                mycount_comp[0] ++;
            } else if(k<Numbers[0]+Numbers[1]){
                mycount_comp[1] ++;
            } else if(k<Numbers[0]+Numbers[1]+Numbers[2]){
                mycount_comp[2] ++;
            } else {
                mycount_comp[3] ++;
            }
            mycount ++;
        }
    }
    // dprintlmpi(mycount);
    // dprintlmpi(Pall.Ntotal_t);
    // dprintlmpi(mycount_comp[0]);
    // exit(0);
#endif

    Pall.Ntotal = mycount;
    Pall.Nhydro = mycount_comp[0];
    Pall.NDM = mycount_comp[1];
    Pall.Nstars = mycount_comp[2]+mycount_comp[3];

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);
    GenerateStructPstar(Pall.Nstars);

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
        Pbody[i]->Type = TypeHydro;
    }
    for(int i=0;i<Pall.Nstars;i++){
        Pbody[i+Pall.Nhydro]->Baryon = (void *)(Pstar[i]);
        Pstar[i]->Body = Pbody[i+Pall.Nhydro];
        Pbody[i+Pall.Nhydro]->Type = TypeStar;
    }

    // dprintlmpi(mycount);
    // dprintlmpi(mycount_comp[0]);
    // dprintlmpi(mycount_comp[1]);
    // dprintlmpi(mycount_comp[2]);
    // dprintlmpi(mycount_comp[3]);

#ifdef FOR_DWARF //{
    double Eps[] = {2.5*PC_CGS/Pall.UnitLength,5.0*PC_CGS/Pall.UnitLength,2.5*PC_CGS/Pall.UnitLength,2.5*PC_CGS/Pall.UnitLength};
    for(int i=0;i<4;i++)
        Eps[i] *= cbrt(10.0);
#else // FOR_DWARF //}//{
    double Eps[] = {25.0*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength};
#endif // FOR_DWARF //}

    //double4 Eps[] = {1.25*PC_CGS/Pall.UnitLength,2.5*PC_CGS/Pall.UnitLength,1.25*PC_CGS/Pall.UnitLength,1.25*PC_CGS/Pall.UnitLength};
    //double Eps[] = {2.5*PC_CGS/Pall.UnitLength,5.0*PC_CGS/Pall.UnitLength,2.5*PC_CGS/Pall.UnitLength,2.5*PC_CGS/Pall.UnitLength};
    //double Eps[] = {100.0*PC_CGS/Pall.UnitLength,200.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {80.0*PC_CGS/Pall.UnitLength,160.0*PC_CGS/Pall.UnitLength,80.0*PC_CGS/Pall.UnitLength,80.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {25.0*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {50.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {100.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength,100.0*PC_CGS/Pall.UnitLength};
    for(int i=0;i<NProcs;i++){
        if(MyID == i){
            ReadDiceLine(fname,Numbers[0]+Numbers[1]+Numbers[2]+Numbers[3],mycount_comp,Eps);
        }
        //MPI_Barrier(MPI_COMM_WORLD);
    }

    ActivateAllparticles();


#if 0
    {

        FILE *fp;
        char fname[MaxCharactersInLine];
        sprintf(fname,"tmpInit.%02d",MyID);
        FileOpen(fp,fname,"w");
        for(int i=0;i<Pall.Nhydro;i++){
            fprintf(fp,"%ld %g %g %g %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                    PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                    PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2],
                    Phydro[i]->Kernel,Phydro[i]->U,PhydroMass(i));
        }
        fclose(fp);
        fflush(NULL);
        MPI_Barrier(MPI_COMM_WORLD);

        if(MPIGetMyID() == MPI_ROOT_RANK){
            system("cat tmpInit.* | sort -n > Init.dat");
            fflush(NULL);
            system("rm -rf ./tmpInit.*");
            fflush(NULL);
        }

    }
#endif



    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.e0;

    Pall.Ns = 128;
    Pall.Npm = 8;

    Pall.hubble = 0.01*70;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

#ifdef FOR_DWARF //{
    Pall.TEnd = 0.4*GIGAYEAR_CGS/Pall.UnitTime;
#else // FOR_DWARF //}//{
    Pall.TEnd = 1*GIGAYEAR_CGS/Pall.UnitTime;
#endif // FOR_DWARF//}
    //Pall.TEnd = 2*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/200.0; 
    //Pall.OutPutInterval = Pall.TEnd/400.0; 

    if(MPIGetMyID() == MPI_ROOT_RANK)
        MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/DICE.ASCII");
    strcpy(Pall.BaseFileName,"./data/DICE");
    strcpy(Pall.RestartFileName,"./data/DICE.dump");

    return;
}


void ReadDiceBinary(char fname[]){

    StructureSizeReport();
#ifdef USE_CELIB
    InitializeStellarFeedback();
#endif

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    
    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MPIGetMyID());

    int Numbers[4] = {0,0,0,0};
    double Units[4] = {0.e0,0.e0,0.e0,0.e0};
    FILE *fp;
    FileOpen(fp,fname,"r");
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fread(Numbers,1,sizeof(int),fp);
        fread(Numbers+1,1,sizeof(int),fp);
        fread(Numbers+2,1,sizeof(int),fp);
        fread(Numbers+3,1,sizeof(int),fp);
        fread(Units,1,sizeof(double),fp);
        fread(Units+1,1,sizeof(double),fp);
        fread(Units+2,1,sizeof(double),fp);
        fread(Units+3,1,sizeof(double),fp);
        /*
        fscanf(fp,"%d %d %d %d",Numbers,Numbers+1,Numbers+2,Numbers+3);
        fscanf(fp,"%le %le %le %le",Units,Units+1,Units+2,Units+3);
        */
    }
    fclose(fp);
    MPI_Bcast(&Numbers,4,MPI_INT,MPI_ROOT_RANK,MPI_COMM_WORLD);
    MPI_Bcast(&Units,4,MPI_DOUBLE,MPI_ROOT_RANK,MPI_COMM_WORLD);

    // Set Shuffule/Lock information
    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"Particle numbers, Ngas,Ndisk,Nhalo,Nspheroid = %d %d %d %d\n",
                Numbers[0],Numbers[1],Numbers[2],Numbers[3]);
        fprintf(stderr,"Simulation units, Mass,Length,Time = %g [g],  %g [cm], %g [s]\n",
                Units[0],Units[1],Units[2]);
    }
    Pall.Ntotal_t = Numbers[0]+Numbers[1]+Numbers[2]+Numbers[3];


    Pall.UnitMass = Units[0];
    Pall.UnitLength = Units[1];
    Pall.UnitTime = Units[2];
    Pall.TCMB = CMB_TEMPERATURE;

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    //SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
    SetVariableViscosityParameters(VARIABLE_VISCOSITY_MIN,VARIABLE_VISCOSITY_MAX,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    int mycount_comp[4] = {0,0,0,0};
    int mycount = 0;

    for(int k=0;k<Pall.Ntotal_t;k++){
        if(k%NProcs == MyID){
            if(k<Numbers[0]){
                mycount_comp[0] ++;
            } else if(k<Numbers[0]+Numbers[1]){
                mycount_comp[1] ++;
            } else if(k<Numbers[0]+Numbers[1]+Numbers[2]){
                mycount_comp[2] ++;
            } else {
                mycount_comp[3] ++;
            }
            mycount ++;
        }
    }

    Pall.Ntotal = mycount;
    Pall.Nhydro = mycount_comp[0];
    Pall.NDM = mycount_comp[1];
    Pall.Nstars = mycount_comp[2]+mycount_comp[3];

    int AllocationSize = mycount;
    if(NProcs > 1)
        AllocationSize = (int)(ForAngelsShare*AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);
    GenerateStructPstar(Pall.Nstars);

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
        Pbody[i]->Type = TypeHydro;
    }
    for(int i=0;i<Pall.Nstars;i++){
        Pbody[i+Pall.Nhydro]->Baryon = (void *)(Pstar[i]);
        Pstar[i]->Body = Pbody[i+Pall.Nhydro];
        Pbody[i+Pall.Nhydro]->Type = TypeStar;
    }

    //double Eps[] = {2.5*PC_CGS/Pall.UnitLength,5.0*PC_CGS/Pall.UnitLength,2.5*PC_CGS/Pall.UnitLength,2.5*PC_CGS/Pall.UnitLength};
    //double Eps[] = {2.5*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {25.0*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {5.0*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength};
    double Eps[] = {10.0*PC_CGS/Pall.UnitLength,10.0*PC_CGS/Pall.UnitLength,10.0*PC_CGS/Pall.UnitLength,10.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {10.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength,10.0*PC_CGS/Pall.UnitLength,10.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {10.0*PC_CGS/Pall.UnitLength,50.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength,25.0*PC_CGS/Pall.UnitLength};
    //double Eps[] = {80.0*PC_CGS/Pall.UnitLength,160.0*PC_CGS/Pall.UnitLength,80.0*PC_CGS/Pall.UnitLength,80.0*PC_CGS/Pall.UnitLength};

    DiceRunParameters.SimulationModel = 0;
    DiceRunParameters.MetallicityType = 3;
    DiceRunParameters.UseEstimateXFe = 1;
    if(DiceRunParameters.MetallicityType!=3) // reset DiceRunParameters.UseEstimateXFe
        DiceRunParameters.UseEstimateXFe = 0;
    DiceRunParameters.Eps[0] = Eps[0];
    DiceRunParameters.Eps[1] = Eps[1];
    DiceRunParameters.Eps[2] = Eps[2];
    DiceRunParameters.Eps[3] = Eps[3];

    // Reset parameters.
    LoadDiceRunParameters();


    for(int i=0;i<NProcs;i++){
        if(MyID == i){
            ReadDiceLineBinary(fname,Numbers[0]+Numbers[1]+Numbers[2]+Numbers[3],mycount_comp,DiceRunParameters.Eps);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.e0;

    Pall.Ns = 128;
    Pall.Npm = 8;

    Pall.hubble = 0.01*70;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

    //Pall.TEnd = 1*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TEnd = 2*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID()==MPI_ROOT_RANK)
        fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = Pall.TEnd/200.0; 
    Pall.OutPutInterval = Pall.TEnd/400.0; 

    if(MPIGetMyID() == MPI_ROOT_RANK)
        MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/DICE.ASCII");
    strcpy(Pall.BaseFileName,"./data/DICE");
    strcpy(Pall.RestartFileName,"./data/DICE.dump");

    return;
}


#endif // TASK_DICE_RUN  //}


#ifdef TASK_MW
double CircularVelocityMilkyWayDisk(const double r);
double CircularVelocityMilkyWayHalo(const double r);
double CircularVelocityMilkyWayExpDisk(const double r);
void InitializeParticleDistributionForMilkyWay(const int Number, const double MassFactor){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = KPC_CGS;
    Pall.UnitTime = KPC_CGS/VELOCITY_KMS_CGS;
    Pall.UnitMass = 1.e+10*MSUN_CGS;
    //Pall.UnitLength = MPC_CGS;
    //Pall.UnitTime = 10.e0*GIGAYEAR_CGS;
    //Pall.UnitMass = 1.e+11*MSUN_CGS;

    Pall.TCMB = CMB_TEMPERATURE;
    Pall.hubble = 0.7;
    Pall.Hubble = Pall.hubble*100*(Pall.UnitLength/MPC_CGS)*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);
    //Pall.OmegaM = 1.0;
    //Pall.OmegaL = 0.0; 

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 10.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(VARIABLE_VISCOSITY_MIN,VARIABLE_VISCOSITY_MAX,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    int mycount_body = 0;
    for(int k=0;k<Number;k+=NProcs)
        mycount_body ++;

    Pall.Ntotal = Pall.Nhydro = mycount_body;
    Pall.Ntotal_t = Pall.Nhydro_t = Number;

    PbodySize = Pall.Ntotal;
    PhydroSize = Pall.Nhydro;
    PbodyElementsSize = Pall.Ntotal;
    PhydroElementsSize = Pall.Nhydro;

    dlprintlmpi(Pall.NDM);
    dlprintlmpi(Pall.Ntotal);
    dlprintlmpi(Pall.Nhydro);

    GenerateStructPbody(Pall.Ntotal);
    GenerateStructPhydro(Pall.Nhydro);

    for(int i=0;i<PhydroSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

#define MgasFrac 0.1
    double Mgas = MgasFrac*4.e+10*MSUN_CGS/Pall.UnitMass;
    double Rdisk = 3.5*KPC_CGS/Pall.UnitLength;
    double z0 = 400.0*PC_CGS/Pall.UnitLength;

//#define DiskEdge (3.0)
#define DiskEdge (1.5)
//#define DiskEdge (1.0/3.5)

    double fact = (1-exp(-DiskEdge)*(1+DiskEdge));
    double mass = fact*Mgas/(double)Number;
    gprintlmpi(fact);
    gprintlmpi(MassFactor*mass);

    const static double OmegaM = 0.3;
    const static double delta_th = 340.e0;

    double Mvir = 1.e+12*MSUN_CGS/Pall.UnitMass;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double Rvir = cbrt(Mvir*3.e0/(4.e0*M_PI*RhoCrit*OmegaM*delta_th));
    fprintf(stderr,"Mvir : %g Msun\n",Mvir*Pall.UnitMass/MSUN_CGS);
    fprintf(stderr,"RhoCrit : %g g/cm^3\n",RhoCrit*Pall.UnitMass/CUBE(Pall.UnitLength));

    int EffectiveNumber = (int)((1.e+11*MSUN_CGS/Pall.UnitMass)/mass);
    double epsilon = 4.0*Rvir/sqrt((double)EffectiveNumber);
    epsilon *= (2/2.8);
    epsilon *= pow(0.1, 1./3.);

    // epsilon is set to the Jean length at SF density.
    double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(1.e+3/Pall.ConvertNumberDensityToCGS/0.76));
    //double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(10*SFDensityCrit/Pall.ConvertNumberDensityToCGS/0.76));
    //double epsilon2 = PC_CGS/Pall.UnitLength;
    //double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(SFDensityCrit/Pall.ConvertNumberDensityToCGS/0.76));

    fprintf(stderr,"Number : %d, EffectiveNumber : %d\n",Number,EffectiveNumber);
    fprintf(stderr,"Rvir : %g kpc\n",Rvir*Pall.UnitLength/KPC_CGS);
    fprintf(stderr,"epsilon : %g kpc, mass : %g Msun\n",
            epsilon*Pall.UnitLength/KPC_CGS,mass*Pall.UnitMass/MSUN_CGS);
    fprintf(stderr,"epsilon2 : %g kpc, SF D %g %g\n",
            epsilon2*Pall.UnitLength/KPC_CGS,SFDensityCrit,SFDensityCrit/Pall.ConvertNumberDensityToCGS);

    fprintf(stderr,"epsilon : %g epsilon2 : %g \n",epsilon,epsilon2);
    epsilon = epsilon2;

    fprintf(stderr,"epsilon x 0.1 !\n");
    //epsilon *= 0.1;
    epsilon = 10*PC_CGS/Pall.UnitLength;
    //epsilon = 1*PC_CGS/Pall.UnitLength;

    //exit(0);

    double InitialKernelSize = cbrt((SQ(DiskEdge*Rdisk)*z0)/(double)Number);
    fprintf(stderr,"Initial Kernel Size : %g kpc\n",InitialKernelSize*Pall.UnitLength/KPC_CGS);

    double Uinit = 10000*Pall.ConvertTtoU;

    int Index = MPIGetMyID();
    //for(int k=MPIGetMyID();k<Pall.Ntotal_t;k+=MPIGetNumProcs()){
    for(int i=0;i<Pall.Ntotal;i++){
        double M = mass*(Index+1);
        double F = M/Mgas;
        
        // solve about r.
        double r = ReturnDistanceForExponent(F,2.0*Rdisk);
        double theta = 2.0*M_PI*gsl_rng_uniform(RandomGenerator);

        double Pos[3];
        Pos[0] = r*cos(theta);
        Pos[1] = r*sin(theta);
        Pos[2] = gsl_ran_gaussian(RandomGenerator,z0);


        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = (i*NProcs+MyID);

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        //Pbody[i]->Pos[2] = Pos[2];
        Pbody[i]->Pos[2] = 0.e0;

        Pbody[i]->PosP[0] = Pbody[i]->Pos[0];
        Pbody[i]->PosP[1] = Pbody[i]->Pos[1];
        Pbody[i]->PosP[2] = Pbody[i]->Pos[2];

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = MassFactor*mass;
        Pbody[i]->Eps = epsilon;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = InitialKernelSize;
        //PbodyHydro(i)->Kernel = epsilon;
        PbodyHydro(i)->U = Uinit;

#if (UseSFModelSpawn) 
        PbodyHydro(i)->SpawnMass = Pbody[i]->Mass/(double)MaxSpawnTimes;
#endif

#ifdef USE_CELIB //{
        // CELibSetPrimordialMetallicity(Pbody[i]->Mass,PbodyHydro(i)->Elements);
        //CELibSetSolarMetallicity(Pbody[i]->Mass,PbodyHydro(i)->Elements);
         CELibSetMetallicityWithSolarAbundancePattern(Pbody[i]->Mass,PbodyHydro(i)->Elements,
                0.5*CELibGetMetalFractionForSolarChemicalComposision());
        // CELibSetMetallicityWithSolarAbundancePattern(Pbody[i]->Mass,PbodyHydro(i)->Elements,0.02);
        double MassLightElements = Phydro[i]->Elements[CELibYield_H]+Phydro[i]->Elements[CELibYield_He];
        double Z = (Pbody[i]->Mass-MassLightElements)/Pbody[i]->Mass;
        PbodyHydro(i)->Z = PbodyHydro(i)->ZII = PbodyHydro(i)->ZIa = Z;
#else 
        PbodyHydro(i)->Z   = 0.02;
        PbodyHydro(i)->ZII = 0.02;
        PbodyHydro(i)->ZIa = 0.00;
#endif // USE_CELIB //}

        Index += MPIGetNumProcs();
    }

    ClearGravitationalForce();
    MilkyWayPotentialHaloDisk();

    // add self-gravity
    Index = MPIGetMyID();
    for(int i=0;i<Pall.Ntotal;i++){ 
        double M = MassFactor*mass*(Index+1);

        double r = sqrt(SQ(Pbody[i]->PosP[0])+SQ(Pbody[i]->PosP[1]));
        if(r>0.e0){
            double gravfact = Pall.GravConst*M/CUBE(r);

            Pbody[i]->Acc[0] -= gravfact*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] -= gravfact*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] -= gravfact*Pbody[i]->PosP[2];
        }
        Index += MPIGetNumProcs();
    }

    double cs = sqrt(5./3.*(5./3.-1.0)*Uinit);
    for(int i=0;i<Pall.Ntotal;i++){
        double Acc = sqrt(SQ(Pbody[i]->Acc[0])+SQ(Pbody[i]->Acc[1]));
        double r = sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1]));
        if(r>TINY){
            Pbody[i]->Vel[0] = -sqrt(Acc/r)*Pbody[i]->Pos[1];
            Pbody[i]->Vel[1] = +sqrt(Acc/r)*Pbody[i]->Pos[0];
        } else {
            Pbody[i]->Vel[0] = gsl_ran_gaussian(RandomGenerator,TINY);
            Pbody[i]->Vel[1] = gsl_ran_gaussian(RandomGenerator,TINY);
        }
        Pbody[i]->Vel[2] = gsl_ran_gaussian(RandomGenerator,cs);
        Pbody[i]->Pos[2] = gsl_ran_gaussian(RandomGenerator,z0);
        Pbody[i]->PosP[2] = Pbody[i]->PosP[2];
    }

    ActivateAllparticles();

#if 0
    {
    FILE *fp,*fpr;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"MW.EXP.Init.%02d.%02d",NProcs,MyID);
    FileOpen(fp ,fname,"w");
    sprintf(fname,"MW.EXP.Init.R.%02d.%02d",NProcs,MyID);
    FileOpen(fpr ,fname,"w");
    double M = 0.0; 
    for(int i=0;i<Pall.Nhydro;i++){
        /*
                */
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                //Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],
                Pall.ConvertUtoT*PbodyHydro(i)->U,Pbody[i]->Mass);
        M += PhydroMass(i);
        fprintf(fpr,"%g %g\n",sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1])),M);
    }
    fclose(fp);
    fclose(fpr);
    //exit(0);
    }
#endif

#if 0
    {
    FILE *fp;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"MW.EXP.CircularProfile.dat");
    FileOpen(fp ,fname,"w");
    double dr = 0.1*Rvir*KPC_CGS/Pall.UnitLength/200.0;
    double r = dr;
    while(r<0.1*Rvir*KPC_CGS/Pall.UnitLength){
        fprintf(fp,"%g %g %g %g %g %g\n",r,CircularVelocityMilkyWayHaloDisk(r),
                CircularVelocityMilkyWayHalo(r),CircularVelocityMilkyWayDisk(r),
                CircularVelocityMilkyWayExpDisk(r),
                CircularVelocityMilkyWayHaloDisk(r)*(Pall.UnitLength/Pall.UnitTime)/1.e+5);
        r += dr;
    }
    fclose(fp);
    exit(0);
    }
#endif

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    // hydro parameters
    // Pall.Ns = 32;
    // Pall.Npm = 2;
    Pall.Ns = 128;
    Pall.Npm = 8;


    // hydro parameters

    Pall.TEnd = 1*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MyID == MPI_ROOT_RANK)
        fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
            Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.001*GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/MW.ASCII");
    strcpy(Pall.BaseFileName,"./data/MW");
    strcpy(Pall.RestartFileName,"./data/MW.dump");



#undef DiskEdge
#undef MgasFrac
    return ;
}

void InitializeParticleDistributionForMilkyWayWithInvRProf(const int Number, const double MassFactor){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = KPC_CGS;
    Pall.UnitTime = KPC_CGS/VELOCITY_KMS_CGS;
    Pall.UnitMass = 1.e+10*MSUN_CGS;

    Pall.TCMB = CMB_TEMPERATURE;
    Pall.hubble = 0.7;
    Pall.Hubble = Pall.hubble*100*(Pall.UnitLength/MPC_CGS)*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 10.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    int mycount_body = 0;
    for(int k=0;k<Number;k+=NProcs)
        mycount_body ++;

    Pall.Ntotal = Pall.Nhydro = mycount_body;
    Pall.Ntotal_t = Pall.Nhydro_t = Number;

    PbodySize = Pall.Ntotal;
    PhydroSize = Pall.Nhydro;
    PbodyElementsSize = Pall.Ntotal;
    PhydroElementsSize = Pall.Nhydro;

    dlprintlmpi(Pall.NDM);
    dlprintlmpi(Pall.Ntotal);
    dlprintlmpi(Pall.Nhydro);

    GenerateStructPbody(Pall.Ntotal);
    GenerateStructPhydro(Pall.Nhydro);

    for(int i=0;i<PhydroSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double Mgas = 3.5e+9*MSUN_CGS/Pall.UnitMass;
    double Rdisk = 2*3.5*KPC_CGS/Pall.UnitLength;
    double z0 = 400.0*PC_CGS/Pall.UnitLength;

#define DiskEdge (1.5) // in unit of the disk scale length.
    double fact = (DiskEdge*7.0)/17.0;  // 17 kpc: Size of HI disk, Nakanishi Sofue 2003.
    double mass = fact*Mgas/(double)Number;
    gprintlmpi(fact);
    gprintlmpi(mass*Pall.UnitMass/MSUN_CGS);
    gprintlmpi(MassFactor*mass);

    const static double OmegaM = 0.3;
    const static double delta_th = 340.e0;

    double Mvir = 1.e+12*MSUN_CGS/Pall.UnitMass;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double Rvir = cbrt(Mvir*3.e0/(4.e0*M_PI*RhoCrit*OmegaM*delta_th));
    fprintf(stderr,"Mvir : %g Msun\n",Mvir*Pall.UnitMass/MSUN_CGS);
    fprintf(stderr,"RhoCrit : %g g/cm^3\n",RhoCrit*Pall.UnitMass/CUBE(Pall.UnitLength));

    int EffectiveNumber = (int)((1.e+11*MSUN_CGS/Pall.UnitMass)/mass);
    double epsilon = 4.0*Rvir/sqrt((double)EffectiveNumber);
    epsilon *= (2/2.8);
    epsilon *= pow(0.1, 1./3.);


    // epsilon is set to the Jean length at SF density.
    double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(1.e+3/Pall.ConvertNumberDensityToCGS/0.76));
    //double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(10*SFDensityCrit/Pall.ConvertNumberDensityToCGS/0.76));
    //double epsilon2 = PC_CGS/Pall.UnitLength;
    //double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(SFDensityCrit/Pall.ConvertNumberDensityToCGS/0.76));

    fprintf(stderr,"Number : %d, EffectiveNumber : %d\n",Number,EffectiveNumber);
    fprintf(stderr,"Rvir : %g kpc\n",Rvir*Pall.UnitLength/KPC_CGS);
    fprintf(stderr,"epsilon : %g kpc, mass : %g Msun\n",
            epsilon*Pall.UnitLength/KPC_CGS,mass*Pall.UnitMass/MSUN_CGS);
    fprintf(stderr,"epsilon2 : %g kpc, SF D %g %g\n",
            epsilon2*Pall.UnitLength/KPC_CGS,SFDensityCrit,SFDensityCrit/Pall.ConvertNumberDensityToCGS);

    fprintf(stderr,"epsilon : %g epsilon2 : %g \n",epsilon,epsilon2);
    epsilon = epsilon2;


    //fprintf(stderr,"epsilon x 0.1 !\n");
    //epsilon *= 0.1;
    fprintf(stderr,"Use Const Eps = 10 pc\n");
    epsilon = 10.0*PC_CGS/Pall.UnitLength;

    //exit(0);

    double InitialKernelSize = cbrt((SQ(DiskEdge*Rdisk)*z0)/(double)Number);
    fprintf(stderr,"Initial Kernel Size : %g pc\n",InitialKernelSize*Pall.UnitLength/PC_CGS);

    double Uinit = 10000*Pall.ConvertTtoU;

    int Index = MPIGetMyID();
    //for(int k=MPIGetMyID();k<Pall.Ntotal_t;k+=MPIGetNumProcs()){
    for(int i=0;i<Pall.Ntotal;i++){
        double M = mass*(Index+1);
        double F = M/(fact*Mgas);
        
        // solve about r.
        double r = ReturnDistanceForInvR(F,1.5*Rdisk);
        double theta = 2.0*M_PI*gsl_rng_uniform(RandomGenerator);

        double Pos[3];
        Pos[0] = r*cos(theta);
        Pos[1] = r*sin(theta);
        Pos[2] = gsl_ran_gaussian(RandomGenerator,z0);


        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = (i*NProcs+MyID);

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        //Pbody[i]->Pos[2] = Pos[2];
        Pbody[i]->Pos[2] = 0.e0;

        Pbody[i]->PosP[0] = Pbody[i]->Pos[0];
        Pbody[i]->PosP[1] = Pbody[i]->Pos[1];
        Pbody[i]->PosP[2] = Pbody[i]->Pos[2];

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = MassFactor*mass;
        Pbody[i]->Eps = epsilon;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = InitialKernelSize;
        //PbodyHydro(i)->Kernel = epsilon;
        PbodyHydro(i)->U = Uinit;

#if (UseSFModelSpawn) 
        PbodyHydro(i)->SpawnMass = Pbody[i]->Mass/(double)MaxSpawnTimes;
#endif
        PbodyHydro(i)->Z   = 0.02;
        PbodyHydro(i)->ZII = 0.02;
        PbodyHydro(i)->ZIa = 0.00;
        Index += MPIGetNumProcs();
    }
    fprintf(stderr,"Particle Mass = %g Msun.\n",Pbody[0]->Mass*Pall.UnitMass/MSUN_CGS);

    ClearGravitationalForce();
    MilkyWayPotentialHaloDisk();

    // add self-gravity
    Index = MPIGetMyID();
    for(int i=0;i<Pall.Ntotal;i++){ 
        double M = MassFactor*mass*(Index+1);

        double r = sqrt(SQ(Pbody[i]->PosP[0])+SQ(Pbody[i]->PosP[1]));
        if(r>0.e0){
            double gravfact = Pall.GravConst*M/CUBE(r);

            Pbody[i]->Acc[0] -= gravfact*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] -= gravfact*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] -= gravfact*Pbody[i]->PosP[2];
        }
        Index += MPIGetNumProcs();
    }

    double cs = sqrt(5./3.*(5./3.-1.0)*Uinit);
    for(int i=0;i<Pall.Ntotal;i++){
        double Acc = sqrt(SQ(Pbody[i]->Acc[0])+SQ(Pbody[i]->Acc[1]));
        double r = sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1]));
        if(r>TINY){
            Pbody[i]->Vel[0] = -sqrt(Acc/r)*Pbody[i]->Pos[1];
            Pbody[i]->Vel[1] = +sqrt(Acc/r)*Pbody[i]->Pos[0];
        } else {
            Pbody[i]->Vel[0] = gsl_ran_gaussian(RandomGenerator,TINY);
            Pbody[i]->Vel[1] = gsl_ran_gaussian(RandomGenerator,TINY);
        }
        Pbody[i]->Vel[2] = gsl_ran_gaussian(RandomGenerator,cs);
        Pbody[i]->Pos[2] = gsl_ran_gaussian(RandomGenerator,z0);
        Pbody[i]->PosP[2] = Pbody[i]->PosP[2];
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

#if 0
    {
    FILE *fp,*fpr;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"Init.%02d.%02d",NProcs,MyID);
    FileOpen(fp ,fname,"w");
    sprintf(fname,"Init.R.%02d.%02d",NProcs,MyID);
    FileOpen(fpr ,fname,"w");
    double M = 0.0; 
    for(int i=0;i<Pall.Nhydro;i++){
        /*
                */
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                //Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],
                Pall.ConvertUtoT*PbodyHydro(i)->U,Pbody[i]->Mass);
        M += Pbody[i]->Mass;
        fprintf(fpr,"%g %g\n",sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1])),M);
    }
    fclose(fp);
    fclose(fpr);
    //exit(0);
    }
#endif

#if 0
    {
    FILE *fp;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"CircularProfile.dat");
    FileOpen(fp ,fname,"w");
    double dr = 0.1*Rvir*KPC_CGS/Pall.UnitLength/200.0;
    double r = dr;
    while(r<0.1*Rvir*KPC_CGS/Pall.UnitLength){
        fprintf(fp,"%g %g %g %g %g %g\n",r,CircularVelocityMilkyWayHaloDisk(r),
                CircularVelocityMilkyWayHalo(r),CircularVelocityMilkyWayDisk(r),
                CircularVelocityMilkyWayExpDisk(r),
                CircularVelocityMilkyWayHaloDisk(r)*(Pall.UnitLength/Pall.UnitTime)/1.e+5);
        r += dr;
    }
    fclose(fp);
    exit(0);
    }
#endif

    Pall.RunStatus = NewSimulation;

    // hydro parameters
    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}
    // hydro parameters

    Pall.TEnd = 1*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.AdaptiveSofteningFactor = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.01*GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/MW.ASCII");
    strcpy(Pall.BaseFileName,"./data/MW");
    strcpy(Pall.RestartFileName,"./data/MW.dump");

#undef DiskEdge 
    return ;
}

static double ReturnDistanceForNFW(const double F, const double Cnfw){

#define Rstart (0)
#define Rend   (4.0)

    double Diff;

    double rstart = Rstart;
    double rend = Rend;
    double fstart =  1.e0/(1.e0 + Cnfw*Rstart) - 1.e0 + log(1.e0 + Cnfw*Rstart) - F;
    double fend = 1.e0/(1.e0 + Cnfw*Rend) - 1.e0 + log(1.e0 + Cnfw*Rend)- F;

    double r = 0.5*(rend+rstart);
    double f = 1.e0/(1.e0 + Cnfw*r) - 1.e0 + log(1.e0 + Cnfw*r) - F;

    int Iteration = 0;
    do{
        f = 1.e0/(1.e0 + Cnfw*r) - 1.e0 + log(1.e0 + Cnfw*r) - F;

        if(f*fstart < 0.0){
            fend = f;
            rend = r;
            r = 0.5*(rstart+rend);
            Diff = fabs(rend-rstart);
        } else {
            fstart = f;
            rstart = r;
            r = 0.5*(rstart+rend);
            Diff = fabs(rend-rstart);
        }
        if(Iteration > MaxIteration)
            exit(1);

        Iteration ++;
    } while ( Diff > 1.e-6*Rend );

    return r;
}

static double ReturnMilkyWayVirialTemperature(void){

    /* Parameters for NFW halo */
    const static double Cnfw = 12.e0;
    const static double OmegaM = 0.3;
    const static double delta_th = 340.e0;

    double Mvir = 1.e+12*MSUN_CGS/Pall.UnitMass;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double Rvir = cbrt(Mvir*3.e0/(4.e0*M_PI*RhoCrit*OmegaM*delta_th));
    double rs = Rvir/Cnfw;
    /* Parameters for NFW halo */

    double delta_0 = (delta_th*OmegaM/3.e0)*(CUBE(Cnfw)/(log(1.e0+Cnfw)-Cnfw/(1.e0+Cnfw)));
    double rho0 = RhoCrit*delta_0;
    double Etot = -16.e0*Pall.GravConst*SQ(M_PI)*SQ(rho0)*pow(rs,5.e0)
        *(Cnfw*(2.e0+Cnfw)-2.e0*(1.e0+Cnfw)*log(1.e0+Cnfw))/(2.e0*SQ(1.e0+Cnfw));
    double Tvir = -Pall.ConvertUtoT*(0.5*Etot/Mvir);

    return Tvir;
}

void InitializeParticleDistributionForMilkyWayWithHaloGas(const int Number){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.UnitLength = KPC_CGS;
    Pall.UnitTime = KPC_CGS/VELOCITY_KMS_CGS;
    Pall.UnitMass = 1.e+10*MSUN_CGS;

    Pall.TCMB = CMB_TEMPERATURE;
    Pall.hubble = 0.7;
    Pall.Hubble = Pall.hubble*100*(Pall.UnitLength/MPC_CGS)*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 10.e0;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    int mycount_body = 0;
    for(int k=0;k<Number;k+=NProcs)
        mycount_body ++;

    int Ndisk = mycount_body;
    int Nhalo = mycount_body;

    Pall.Ntotal = Pall.Nhydro = 2*mycount_body;
    Pall.Ntotal_t = Pall.Nhydro_t = 2*Number;

    PbodySize = Pall.Ntotal;
    PhydroSize = Pall.Nhydro;
    PbodyElementsSize = Pall.Ntotal;
    PhydroElementsSize = Pall.Nhydro;

    dlprintlmpi(Pall.NDM);
    dlprintlmpi(Pall.Ntotal);
    dlprintlmpi(Pall.Nhydro);

    GenerateStructPbody(Pall.Ntotal);
    GenerateStructPhydro(Pall.Nhydro);

    for(int i=0;i<PhydroSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

#define MgasFrac 0.1
    double Mgas = MgasFrac*4.e+10*MSUN_CGS/Pall.UnitMass;
    double Rdisk = 3.5*KPC_CGS/Pall.UnitLength;
    double z0 = 400.0*PC_CGS/Pall.UnitLength;

#define DiskEdge (1.5)

    double fact = (1-exp(-DiskEdge)*(1+DiskEdge));
    double mass = fact*Mgas/(double)Number;
    gprintlmpi(fact);
    gprintlmpi(mass);

    const static double OmegaM = 0.3;
    const static double delta_th = 340.e0;

    double Mvir = 1.e+12*MSUN_CGS/Pall.UnitMass;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double Rvir = cbrt(Mvir*3.e0/(4.e0*M_PI*RhoCrit*OmegaM*delta_th));
    fprintf(stderr,"Mvir : %g Msun\n",Mvir*Pall.UnitMass/MSUN_CGS);
    fprintf(stderr,"RhoCrit : %g g/cm^3\n",RhoCrit*Pall.UnitMass/CUBE(Pall.UnitLength));

    int EffectiveNumber = (int)((1.e+11*MSUN_CGS/Pall.UnitMass)/mass);
    double epsilon = 4.0*Rvir/sqrt((double)EffectiveNumber);
    epsilon *= (2/2.8);
    epsilon *= pow(0.1, 1./3.);

    // epsilon is set to the Jean length at SF density.
    //double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(1.e+3/Pall.ConvertNumberDensityToCGS/0.76));
    double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(10*SFDensityCrit/Pall.ConvertNumberDensityToCGS/0.76));
    //double epsilon2 = cbrt((3.0/(4.0*M_PI))*(2.0*32*mass)/(SFDensityCrit/Pall.ConvertNumberDensityToCGS/0.76));

    fprintf(stderr,"Number : %d, EffectiveNumber : %d\n",Number,EffectiveNumber);
    fprintf(stderr,"Rvir : %g kpc\n",Rvir*Pall.UnitLength/KPC_CGS);
    fprintf(stderr,"epsilon : %g kpc, mass : %g Msun\n",
            epsilon*Pall.UnitLength/KPC_CGS,mass*Pall.UnitMass/MSUN_CGS);
    fprintf(stderr,"epsilon2 : %g kpc, SF D %g %g\n",
            epsilon2*Pall.UnitLength/KPC_CGS,SFDensityCrit,SFDensityCrit/Pall.ConvertNumberDensityToCGS);

    fprintf(stderr,"epsilon : %g epsilon2 : %g \n",epsilon,epsilon2);
    epsilon = epsilon2;
    //exit(0);

    double InitialKernelSize = cbrt((SQ(DiskEdge*Rdisk)*z0)/(double)Number);
    fprintf(stderr,"Initial Kernel Size : %g kpc\n",InitialKernelSize*Pall.UnitLength/KPC_CGS);

    double Uinit = 10000*Pall.ConvertTtoU;

    int Index = MPIGetMyID();
    for(int i=0;i<Ndisk;i++){
        double M = mass*(Index+1);
        double F = M/Mgas;
        
        // solve about r.
        double r = ReturnDistanceForExponent(F,2.0*Rdisk);
        double theta = 2.0*M_PI*gsl_rng_uniform(RandomGenerator);

        double Pos[3];
        Pos[0] = r*cos(theta);
        Pos[1] = r*sin(theta);
        Pos[2] = gsl_ran_gaussian(RandomGenerator,z0);


        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = (i*NProcs+MyID);

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        Pbody[i]->Pos[2] = 0.e0;

        Pbody[i]->PosP[0] = Pbody[i]->Pos[0];
        Pbody[i]->PosP[1] = Pbody[i]->Pos[1];
        Pbody[i]->PosP[2] = Pbody[i]->Pos[2];

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = 2.0*mass;
        Pbody[i]->Eps = epsilon;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = InitialKernelSize;
        PbodyHydro(i)->U = Uinit;

#if (UseSFModelSpawn) 
        PbodyHydro(i)->SpawnMass = Pbody[i]->Mass/(double)MaxSpawnTimes;
#endif
        PbodyHydro(i)->Z   = 0.02;
        PbodyHydro(i)->ZII = 0.02;
        PbodyHydro(i)->ZIa = 0.00;
        Index += MPIGetNumProcs();
    }

    // set nfw halo
    double Cnfw = 12.e0;
    double fend = 1/(1 + Cnfw) - 1 + log(1 + Cnfw);
    double dfend = fend/((double)Nhalo);
    Index = MPIGetMyID();
    double Uinit_halo = ReturnMilkyWayVirialTemperature();
    //double Uinit_halo = Pall.ConvertUtoT*Uinit;
    fprintf(stderr,"Virial Temperature = %e \n",Uinit_halo);
    for(int i=Ndisk;i<(Ndisk+Nhalo);i++){
        double F = dfend*(Index+1);
        
        // solve about r.
        double r = Rvir*ReturnDistanceForNFW(F,Cnfw);
        double z = r*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        double R = sqrt(SQ(r)-SQ(z));
        double theta = 2.0*M_PI*gsl_rng_uniform(RandomGenerator);

        double Pos[3];
        Pos[0] = R*cos(theta);
        Pos[1] = R*sin(theta);
        Pos[2] = z;


        Pbody[i]->Active = ON;
        Pbody[i]->Use = ON;
        Pbody[i]->Type = TypeHydro;
        Pbody[i]->GlobalID = (i*NProcs+MyID);

        Pbody[i]->Pos[0] = Pos[0];
        Pbody[i]->Pos[1] = Pos[1];
        Pbody[i]->Pos[2] = Pos[2];

        Pbody[i]->PosP[0] = Pbody[i]->Pos[0];
        Pbody[i]->PosP[1] = Pbody[i]->Pos[1];
        Pbody[i]->PosP[2] = TINY;

        Pbody[i]->Vel[0] = TINY;
        Pbody[i]->Vel[1] = TINY;
        Pbody[i]->Vel[2] = TINY;

        Pbody[i]->Mass = 2.0*mass;
        Pbody[i]->Eps = epsilon;

        PbodyHydro(i)->Use = ON;
        PbodyHydro(i)->Kernel = InitialKernelSize;
        PbodyHydro(i)->U = Pall.ConvertTtoU*Uinit_halo;

#if (UseSFModelSpawn) 
        PbodyHydro(i)->SpawnMass = Pbody[i]->Mass/(double)MaxSpawnTimes;
#endif
        PbodyHydro(i)->Z   = 0.02;
        PbodyHydro(i)->ZII = 0.02;
        PbodyHydro(i)->ZIa = 0.00;
        Index += MPIGetNumProcs();
    }

    ClearGravitationalForce();
    MilkyWayPotentialHaloDisk();

    double cs = sqrt(5./3.*(5./3.-1.0)*Uinit);
    for(int i=0;i<Ndisk;i++){
        double Acc = sqrt(SQ(Pbody[i]->Acc[0])+SQ(Pbody[i]->Acc[1]));
        double r = sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1]));
        if(r>TINY){
            Pbody[i]->Vel[0] = -sqrt(Acc/r)*Pbody[i]->Pos[1];
            Pbody[i]->Vel[1] = +sqrt(Acc/r)*Pbody[i]->Pos[0];
        } else {
            Pbody[i]->Vel[0] = gsl_ran_gaussian(RandomGenerator,TINY);
            Pbody[i]->Vel[1] = gsl_ran_gaussian(RandomGenerator,TINY);
        }
        Pbody[i]->Vel[2] = gsl_ran_gaussian(RandomGenerator,cs);
        Pbody[i]->PosP[2] = Pbody[i]->Pos[2] = gsl_ran_gaussian(RandomGenerator,z0);
    }


    cs = sqrt(5./3.*(5./3.-1.0)*Uinit_halo);
    for(int i=Ndisk;i<(Ndisk+Nhalo);i++){
        double Acc = sqrt(SQ(Pbody[i]->Acc[0])+SQ(Pbody[i]->Acc[1]));
        double r = sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1]));
        if(r>TINY){
            Pbody[i]->Vel[0] = -sqrt(Acc/r)*Pbody[i]->Pos[1];
            Pbody[i]->Vel[1] = +sqrt(Acc/r)*Pbody[i]->Pos[0];
        } else {
            Pbody[i]->Vel[0] = gsl_ran_gaussian(RandomGenerator,TINY);
            Pbody[i]->Vel[1] = gsl_ran_gaussian(RandomGenerator,TINY);
        }
        Pbody[i]->Vel[2] = gsl_ran_gaussian(RandomGenerator,cs);
        Pbody[i]->PosP[2] = Pbody[i]->Pos[2];
    }


#if 1
    FILE *fp;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"Init.%02d.%02d",NProcs,MyID);
    FileOpen(fp ,fname,"w");
    //double M = 0.0; 
    for(int i=0;i<Pall.Nhydro;i++){
        /*
                */
        fprintf(fp,"%ld %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                //Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2],
                Pall.ConvertUtoT*PbodyHydro(i)->U,Pbody[i]->Mass);
        //M += Pbody[i]->Mass;
        //fprintf(fp,"%g %g\n",sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1])),M);
    }
    fclose(fp);
    //exit(0);
#endif

#if 0
    FILE *fp;
    char fname[MaxCharactersInLine]; 
    sprintf(fname,"CircularProfile.dat");
    FileOpen(fp ,fname,"w");
    double dr = 0.1*Rvir*KPC_CGS/Pall.UnitLength/200.0;
    double r = dr;
    while(r<0.1*Rvir*KPC_CGS/Pall.UnitLength){
        fprintf(fp,"%g %g %g %g %g %g\n",r,CircularVelocityMilkyWayHaloDisk(r),
                CircularVelocityMilkyWayHalo(r),CircularVelocityMilkyWayDisk(r),
                CircularVelocityMilkyWayExpDisk(r),
                CircularVelocityMilkyWayHaloDisk(r)*(Pall.UnitLength/Pall.UnitTime)/1.e+5);
        r += dr;
    }
    fclose(fp);
    exit(0);
#endif

    Pall.RunStatus = NewSimulation;

    // hydro parameters
    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}
    // hydro parameters

    Pall.TEnd = 1.0*GIGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.AdaptiveSofteningFactor = 1.e0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.05*GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/MW.ASCII");
    strcpy(Pall.BaseFileName,"./data/MW");
    strcpy(Pall.RestartFileName,"./data/MW.dump");

    return ;
}

static inline double __attribute__((always_inline)) fx(const double x){
    return (log(1.0+x) - x/(1.e0+x));
}

void ExternalPotentialForNFW(const double M200){

    //double M200 = pow(10.0,(double)(M200_Log10M-10));
    //dprintlmpi(M200_Log10M);
    //gprintlmpi(M200);
    //exit(0);
    double Cnfw = 10.e0;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double delta_c = 200 * CUBE(Cnfw)/(3*fx(Cnfw));
    double rho0 = RhoCrit * delta_c;
    double R200 = pow(M200*3/(200*RhoCrit*4*M_PI), 1.e0/3.e0);
    double rs = R200/Cnfw;
    double Softening = 10*PC_CGS/Pall.UnitLength;
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            double r = sqrt(NORM2(Pbody[i]->PosP)+SQ(0.1*rs))+Softening;
            // double r = sqrt(NORM2(Pbody[i]->PosP)+SQ(0.1*rs));
            double x = r/rs;
            double Mass = 4*M_PI*rho0*CUBE(rs)*fx(x);
            double gravfact = Pall.GravConst*Mass/CUBE(r);

            Pbody[i]->Acc[0] -= gravfact*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] -= gravfact*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] -= gravfact*Pbody[i]->PosP[2];
        }
    }

    return ;
}

#define MiyamotoNagaiDiskScaleLength (3.5) // in kpc
#define MiyamotoNagaiDiskScaleHeight (400.0) // in pc
void MilkyWayPotentialHaloDisk(void){

    /* Parameters for NFW halo */
    const static double Cnfw = 12.e0;
    const static double OmegaM = 0.3;
    const static double delta_th = 340.e0;

    double Mvir = 1.e+12*MSUN_CGS/Pall.UnitMass;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double Rvir = cbrt(Mvir*3.e0/(4.e0*M_PI*RhoCrit*OmegaM*delta_th));
    double rs = Rvir/Cnfw;
    /* Parameters for NFW halo */

    /* Parameters for Miyamoto-Nagai disk */
    double Mdisk = 4.e+10*MSUN_CGS/Pall.UnitMass;
    double R0 = MiyamotoNagaiDiskScaleLength*KPC_CGS/Pall.UnitLength;
    double z0 = MiyamotoNagaiDiskScaleHeight*PC_CGS/Pall.UnitLength;
    /* Parameters for Miyamoto-Nagai disk */

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            // NFW halo 
            double r = sqrt(NORM2(Pbody[i]->PosP)+SQ(0.1*rs));
            double x = r/rs;
            double Mass = Mvir*fx(x)/fx(Cnfw);
            double gravfact = Pall.GravConst*Mass/CUBE(r);

            Pbody[i]->Acc[0] -= gravfact*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] -= gravfact*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] -= gravfact*Pbody[i]->PosP[2];

            // Miyamoto Nagai
            double Rdisk = sqrt(SQ(Pbody[i]->PosP[0])+SQ(Pbody[i]->PosP[1]));
            double Zdisk = Pbody[i]->PosP[2];
            double Adisk = R0+sqrt(SQ(Zdisk)+SQ(z0));
            double gravfact_disk = Pall.GravConst*Mdisk/sqrt(CUBE(SQ(Rdisk)+SQ(Adisk)));

            Pbody[i]->Acc[0] -= gravfact_disk*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] -= gravfact_disk*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] -= gravfact_disk*(R0/sqrt(SQ(Zdisk)+SQ(z0))+1.0)*Pbody[i]->PosP[2];
        }
    }

    return ;
}

double CircularVelocityMilkyWayHaloDisk(const double r){

    /* Parameters for NFW halo */
    const static double Cnfw = 12.e0;
    const static double OmegaM = 0.3;
    const static double delta_th = 340.e0;

    double Mvir = 1.e+12*MSUN_CGS/Pall.UnitMass;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double Rvir = cbrt(Mvir*3.e0/(4.e0*M_PI*RhoCrit*OmegaM*delta_th));
    double rs = Rvir/Cnfw;
    /* Parameters for NFW halo */

    /* Parameters for Miyamoto-Nagai disk */
    double Mdisk = 4.e+10*MSUN_CGS/Pall.UnitMass;
    double R0 = MiyamotoNagaiDiskScaleLength*KPC_CGS/Pall.UnitLength;
    double z0 = MiyamotoNagaiDiskScaleHeight*PC_CGS/Pall.UnitLength;
    /* Parameters for Miyamoto-Nagai disk */

    // NFW halo 
    double x = r/rs;
    //double C = Rvir/rs;
    double Mass = Mvir*fx(x)/fx(Cnfw);
    double gravfact_halo = Pall.GravConst*Mass/CUBE(r+0.1*x);
    //double gravfact_halo = Pall.GravConst*Mass/CUBE(r);

    //return (sqrt(Pall.GravConst*Mass/r));

    //fprintf(stderr,"Mvir %g, Mr %g, Rvir %g, rs = %g, x = %g, r = %g \n",
            //Mvir*Pall.UnitMass/MSUN_CGS,Mass*Pall.UnitMass/MSUN_CGS,
            //Rvir*Pall.UnitLength/KPC_CGS,rs*Pall.UnitLength/KPC_CGS,r/rs,r*Pall.UnitLength/KPC_CGS);

    // Miyamoto Nagai
    double Adisk = R0+sqrt(SQ(0.e0)+SQ(z0));
    double gravfact_disk = Pall.GravConst*Mdisk/(CUBE(sqrt((SQ(r)+SQ(Adisk)))));

    // gravfact_disk = 0.e0;
    // gravfact_halo = 0.e0;

    return (sqrt(fabs(gravfact_halo+gravfact_disk))*r);
}

double CircularVelocityMilkyWayHalo(const double r){

    /* Parameters for NFW halo */
    const static double Cnfw = 12.e0;
    const static double OmegaM = 0.3;
    const static double delta_th = 340.e0;

    double Mvir = 1.e+12*MSUN_CGS/Pall.UnitMass;
    double RhoCrit = (3.e0*Pall.Hubble*Pall.Hubble/(8.e0*M_PI*Pall.GravConst));
    double Rvir = cbrt(Mvir*3.e0/(4.e0*M_PI*RhoCrit*OmegaM*delta_th));
    double rs = Rvir/Cnfw;
    /* Parameters for NFW halo */

    // NFW halo 
    double x = r/rs;
    double Mass = Mvir*fx(x)/fx(Cnfw);
    //double gravfact_halo = Pall.GravConst*Mass/CUBE(r+0.1*x);
    double gravfact_halo = Pall.GravConst*Mass/CUBE(r);

    return (sqrt(gravfact_halo)*r);
}

double CircularVelocityMilkyWayDisk(const double r){

    /* Parameters for Miyamoto-Nagai disk */
    double Mdisk = 4.e+10*MSUN_CGS/Pall.UnitMass;
    double r0 = MiyamotoNagaiDiskScaleLength*KPC_CGS/Pall.UnitLength;
    double z0 = 0.5*MiyamotoNagaiDiskScaleHeight*PC_CGS/Pall.UnitLength;
    /* Parameters for Miyamoto-Nagai disk */

    // Miyamoto Nagai
    double Adisk = r0+sqrt(SQ(0.e0)+SQ(z0));
    double gravfact_disk = Pall.GravConst*Mdisk/(CUBE(sqrt((SQ(r)+SQ(Adisk)))));

    return (sqrt(fabs(gravfact_disk))*r);
}

double CircularVelocityMilkyWayExpDisk(const double r){

    /* Parameters for exp disk */
    double Mdisk = 4.e+10*MSUN_CGS/Pall.UnitMass;
    double rd = 3.5*KPC_CGS/Pall.UnitLength;
    double z0 = 400*PC_CGS/Pall.UnitLength;
    /* Parameters for exp disk */
    
    double vc2 = Pall.GravConst*Mdisk*(1.0-(1.0+r/rd)*exp(-r/rd))/r;
    
    return (sqrt(vc2));
}

#if 0
void ParticleRemover(const double Radius){

    // remove very distance particles
    int NLoop = Pall.Ntotal;
    int counter = 0;
    for(int i=0;i<NLoop;i++){
        if(NORM(Pbody[i]->Pos)>Radius){
            Pbody[i]->Use = OFF;
            Pall.Ntotal --;
            if(Pbody[i]->Type == TypeHydro){
                PbodyHydro(i)->Use = OFF;
                Pall.Nhydro --;
            }else if(Pbody[i]->Type == TypeStar){
                PbodyStar(i)->Use = OFF;
                Pall.Nstars --;
            }else if(Pbody[i]->Type == TypeDM){
                Pall.NDM --;
            }
            counter ++;
        }
    }
    //dprintlmpi(counter);

    int sum_counter = 0;
    MPI_Allreduce(&counter,&sum_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(sum_counter > 0){
        ReConnectPointers();
        UpdateTotalNumber();
    }

    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"%d particles are removed in this simulation.\n",sum_counter);
        fprintf(stderr," Removal rule is that R > %g kpc.\n", Radius*Pall.UnitLength/KPC_CGS);
    }

    return ;
}
#endif

void CheckKennicutLaw(const int mode){

    FILE *fp;
    if(mode == 0){ //Initialize
        if(MPIGetMyID()==MPI_ROOT_RANK){
            FileOpen(fp,"./data/Kennicut.dat","w");
            fclose(fp);
        }
    } else { // Check;
#define TSNIIend    (4.5e+7) // yr

        static double Radius[] = {1.0,3.0,5.0,10.0};
        int NBins = sizeof(Radius)/sizeof(double);

        double Mgas[NBins],GlobalMgas[NBins];
        double Mstar[NBins],GlobalMstar[NBins];

        for(int k=0;k<NBins;k++){
            //get gas mass in parallel.
            Mgas[k] = 0.e0;
            for(int i=0;i<Pall.Nhydro;i++){
                double r = Pall.UnitLength*NORM(PhydroPosP(i))/KPC_CGS;
                if(r<Radius[k]){
                    Mgas[k] += PhydroMass(i); 
                }
            }
            //get newly born stellar mass in parallel.
            Mstar[k] = 0.e0;
            for(int i=0;i<Pall.Nstars;i++){
                double r = Pall.UnitLength*NORM(PstarPosP(i))/KPC_CGS;
                if((r<Radius[k])
                    &&(Pall.TCurrent-Pstar[i]->FormationTime < TSNIIend*YEAR_CGS/Pall.UnitTime)){
                    Mstar[k] += PstarMass(i); 
                }
            }
        }
        MPI_Allreduce(Mgas,GlobalMgas,NBins,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(Mstar,GlobalMstar,NBins,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        if(MPIGetMyID()==MPI_ROOT_RANK){
            // write data
            FileOpen(fp,"./data/Kennicut.dat","a");
            // t,z Mg0, Msf0, Mg1, Msf1, Mg2, Msf2, Mg3, Msf3,
            double mgfact = (Pall.UnitMass/MSUN_CGS)/(4*M_PI)*1.e-6; //Msun/pc^2
            double msfact = (Pall.UnitMass/MSUN_CGS)/(4*M_PI)/(TSNIIend); //Msun/kpc-2/yr
            fprintf(fp,"%g %g %g %g %g %g %g %g %g %g\n",Pall.TCurrent,Pall.Redshift,
                mgfact*Mgas[0]/SQ(Radius[0]),msfact*Mstar[0]/SQ(Radius[0]),
                mgfact*Mgas[1]/SQ(Radius[1]),msfact*Mstar[1]/SQ(Radius[1]),
                mgfact*Mgas[2]/SQ(Radius[2]),msfact*Mstar[2]/SQ(Radius[2]),
                mgfact*Mgas[3]/SQ(Radius[3]),msfact*Mstar[3]/SQ(Radius[3]));
            fclose(fp);
        }
    }

    return; 
}
#endif //TASK_MW

#ifdef TASK_AGNTORUS
#define AGN_TORUS_RUN_LIVE_BH ON 
/*
 * This initial condition is based on Wada Norman (2002), Obscuring Material
 * around Seyfert Nuclei with Starburst, ApJ 566, L21-24.
 * The potential is \nabla \Phi_ext + \nabla \Phi_BH.
 * \Phi_ext = (-27/4)^{1/2} v_c^2/(r^2+a^2)^{1/2}, where a = 10pc and
 * vc=100km/s.
 * \Phi_BH = -G M_BH/(r^2+b^2)^{1/2}, where b = 1pc and M_BH = 10^8Msun.
 * The total mass of gas disk is set to be 5x10^7Msun and the scale height is
 * set to be 2.5pc.
 */

extern double Redge_AGN_CGS;
extern double MassBH_AGN_CGS;
extern double MassGas_AGN_CGS;

//#define MGas_AGN    (6.e+6*MSUN_CGS)
#define Vcir_AGN  (147*VELOCITY_KMS_CGS) //cm/s
//#define MBH_AGN (1.3e+7*MSUN_CGS) // Msun
#define Ra_AGN (10*PC_CGS) // R_core for bulge.
#define Ra1_AGN (100*PC_CGS) // R_core for bulge.
#define Ra2_AGN (2500*PC_CGS) // R_core for bulge.
#define Rb_AGN (1*PC_CGS)  // R_core for BH
//#define DISK_EDGE_AGN   (32*PC_CGS) 

//#define VC_ADJUST_AGN (9.e+40)
//#define VC_ADJUST_AGN (2.e+3)
#define VC_ADJUST_AGN (3.098e+3) // PC/YR**2

#define MBULGE_AGN  (1.e+10*MSUN_CGS)
#define RaHernquist_AGN  (450*PC_CGS)

#define USE_HERNQUIST_BULGE_FOR_AGN

static double VcUnit,VcUnit2,FVcUnit2;
static double MBHUnit,GMBHUnit,GMSelfUnit;
static double RaUnit,RaUnit2,RbUnit,RbUnit2; 
static double Ra1Unit,Ra1Unit2,Ra2Unit,Ra2Unit2; 
static double Rmax;
static double GMBulgeUnit,RaHernquistUnit;

void InitializeAGNTorusExternalPotentials(void){

#ifdef USE_HERNQUIST_BULGE_FOR_AGN
    GMBulgeUnit = Pall.GravConst*(MBULGE_AGN/Pall.UnitMass);
    RaHernquistUnit = RaHernquist_AGN/Pall.UnitLength;
#else
    VcUnit  = Vcir_AGN*GetUnitVel();
    VcUnit2 = SQ(VcUnit);
    FVcUnit2 = sqrt(VC_ADJUST_AGN*27.0/4.0)*VcUnit2;
    RaUnit = Ra_AGN/Pall.UnitLength;
    RaUnit2 = SQ(RaUnit);
    RbUnit = Rb_AGN/Pall.UnitLength;
    RbUnit2 = SQ(RbUnit);

    Ra1Unit = Ra1_AGN/Pall.UnitLength;
    Ra1Unit2 = SQ(Ra1Unit);
    Ra2Unit = Ra2_AGN/Pall.UnitLength;
    Ra2Unit2 = SQ(Ra2Unit);

    fprintf(stderr,"--%g %g %g\n",VcUnit,VcUnit2,FVcUnit2);
#endif
    //fprintf(stderr,"--%g %g %g %g %g\n",VcUnit,VcUnit2,FVcUnit2,MBHUnit,GMBHUnit);

    return;
}

double InitializeCircularVelocityFromAGNExternalPotentials(const int index, const double r){

#ifndef USE_HERNQUIST_BULGE_FOR_AGN
    double r2 = SQ(r);
    //double InvRa = 1.e0/sqrt(r2+RaUnit2);
    double InvRa = 1.e0/sqrt(r2+Ra1Unit2)+1.0/sqrt(r2+Ra2Unit2);
    double InvRa3 = CUBE(InvRa);

    return sqrt(FVcUnit2*InvRa3*r2);
#endif
    return 0.0;
}

void AGNForceFromExternalPotentials(void){

    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Active){
            //double InvRa = 1.e0/sqrt(r2+RaUnit2);
#ifdef USE_HERNQUIST_BULGE_FOR_AGN
            double r = NORM(Pbody[i]->PosP);
            double InvRa = 1.0/(r+RaHernquistUnit);
            //double InvRa3 = CUBE(InvRa);
            double InvRa3 = SQ(InvRa)/(r+1.e-10*RaHernquistUnit);
            Pbody[i]->Acc[0] -= (GMBulgeUnit*InvRa3)*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] -= (GMBulgeUnit*InvRa3)*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] -= (GMBulgeUnit*InvRa3)*Pbody[i]->PosP[2];
#else
            double r2 = NORM2(Pbody[i]->PosP);
            double InvRa = 1.e0/sqrt(r2+Ra1Unit2)+1.0/sqrt(r2+Ra2Unit2);
            double InvRa3 = CUBE(InvRa);
            Pbody[i]->Acc[0] -= (FVcUnit2*InvRa3)*Pbody[i]->PosP[0];
            Pbody[i]->Acc[1] -= (FVcUnit2*InvRa3)*Pbody[i]->PosP[1];
            Pbody[i]->Acc[2] -= (FVcUnit2*InvRa3)*Pbody[i]->PosP[2];
#endif
        }
    }
    return;
}

inline static double VcBH_AGN(double R, double EPS){
    //return sqrt(Pall.GravConst*(MBH_AGN/Pall.UnitMass)*SQ(R)*CUBE(1.0/sqrt(SQ(R)+SQ(EPS))));
    return sqrt(Pall.GravConst*(MassBH_AGN_CGS/Pall.UnitMass)*SQ(R)*CUBE(1.0/sqrt(SQ(R)+SQ(EPS))));
}
inline static double VcBulge_AGN(double R){
    return sqrt(Pall.GravConst*(MBULGE_AGN/Pall.UnitMass)*R*(1.0/(SQ(R)+SQ(RaHernquistUnit))));
}
inline static double VcGas_AGN(double R){
    //return sqrt(Pall.GravConst*(MGas_AGN/Pall.UnitMass)*(SQ(R)/SQ(Redge_AGN_CGS/Pall.UnitLength))/R);
    return sqrt(Pall.GravConst*(MassGas_AGN_CGS/Pall.UnitMass)*(SQ(R)/SQ(Redge_AGN_CGS/Pall.UnitLength))/R);
}


void InitAGNTorus(const int Number, const int DistortionFlag){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    /* */
    fprintf(stderr,"%g %g %g\n",
            Redge_AGN_CGS/PC_CGS,
            MassBH_AGN_CGS/MSUN_CGS,
            MassGas_AGN_CGS/MSUN_CGS);

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    //InitializeRandomGenerator(1977+MyID);
    InitializeRandomGenerator(1977);

    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = 1.e+8*MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    InitializeAGNTorusExternalPotentials();

    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pall.Nhydro = Number/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        //Pall.NDM = Number/NProcs;
        //Pall.NDM_t = Pall.NDM*NProcs;
        Pall.Nsink = Pall.Nsink_t = 1;
        Pall.Ntotal = Pall.Nhydro + Pall.Nsink;
        Pall.Ntotal_t = Pall.Nhydro_t + Pall.Nsink_t;
        //Pall.Ntotal = Pall.NDM + Pall.Nsink;
        //Pall.Ntotal_t = Pall.NDM_t + Pall.Nsink_t;
    } else {
        Pall.Nhydro = Number/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        //Pall.NDM = Number/NProcs;
        //Pall.NDM_t = Pall.NDM*NProcs;
        Pall.Nsink = 0;
        Pall.Nsink_t = 1;
        Pall.Ntotal = Pall.Nhydro;
        Pall.Ntotal_t = Pall.Nhydro_t + Pall.Nsink_t;
        //Pall.Ntotal = Pall.NDM;
        //Pall.Ntotal_t = Pall.NDM_t + Pall.Nsink_t;
    }
    InitializeAllActiveParticleNumbers();

    int AllocationSize = Pall.Ntotal; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        GenerateStructPsink(1);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    
    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pbody[AllocationSize-1]->Baryon = (void *)(Psink[0]);
        Psink[0]->Body = Pbody[AllocationSize-1];
    }

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    double mgas = (MassGas_AGN_CGS/Pall.UnitMass)/Pall.Nhydro_t;
    fprintf(stderr,"%g MGas Pall.UntiMass, Mgas = %g [Msun], Nhydro %ld\n",
            mgas,MassGas_AGN_CGS/MSUN_CGS,Pall.Nhydro_t);
    //double mgas = (MGas_AGN/Pall.UnitMass)/Pall.Nhydro_t;
    //fprintf(stderr,"%g MGas Pall.UntiMass, Nhydro %g %g %ld\n",mgas, MGas_AGN/MSUN_CGS,Pall.UnitMass,Pall.Nhydro_t);
    gprintlmpi(mgas);
    //double eps = 1.e0*PC_CGS/Pall.UnitLength;
    double eps = 0.1*PC_CGS/Pall.UnitLength;
    fprintf(stderr,"Eps = %g pc\n",eps*Pall.UnitLength/PC_CGS);
    double Uinit = 1000*Pall.ConvertTtoU;
    //double Uinit = 0.10*Pall.ConvertTtoU;
    double cs = sqrt(Pall.GGm1*Uinit);
    double Zsig = 2.5*PC_CGS/Pall.UnitLength;
    fprintf(stderr,"%g %g %g %g\n",Uinit,cs,Pall.ConvertTtoU,Pall.GGm1);
    //exit(0);
    Rmax = (Redge_AGN_CGS/Pall.UnitLength);
    // Rmax = (DISK_EDGE_AGN/Pall.UnitLength);

#define __AGN_TORUS_TEST__

    int count = 0;
#if 1  // FB test
    int _cccc = 0; int ttt = 0;
    while(_cccc < Pall.Nhydro_t){
        double Pos[3];
        Pos[0] = Rmax*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        Pos[1] = Rmax*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        double r = sqrt(SQ(Pos[0])+SQ(Pos[1]));
        double z = Gaussian();
        double vx = Gaussian();
        double vy = Gaussian();
        double vz = Gaussian();
        ttt ++;
        if((r>1*PC_CGS/Pall.UnitLength) && (r<Rmax)){
            _cccc ++;
            if((_cccc-1)%NProcs != MyID) continue;

            Pbody[count]->Active = ON;
            Pbody[count]->Use = ON;
            Pbody[count]->Type = TypeHydro;
            Pbody[count]->GlobalID = _cccc;

            Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
            Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
            Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Zsig*z;

#if 0
            double Vc = InitializeCircularVelocityFromAGNExternalPotentials(count,r);
            Pbody[count]->Vel[0] = -Vc*Pos[1]/r;
            Pbody[count]->Vel[1] =  Vc*Pos[0]/r;
            Pbody[count]->Vel[2] = 0.e0;

            Pbody[count]->Vel[0] += 0.3*Vc*vx;
            Pbody[count]->Vel[1] += 0.3*Vc*vy;
            Pbody[count]->Vel[2] += 0.3*Vc*vz;
#endif
            Pbody[count]->Vel[0] = Pbody[count]->Vel[1] = Pbody[count]->Vel[2] = 0.0;

            Pbody[count]->Eps = eps;
            Pbody[count]->Mass = mgas;
#if 1
            PbodyHydro(count)->Mass = Pbody[count]->Mass = mgas;
            PbodyHydro(count)->Z = 0.02;
            PbodyHydro(count)->ZII = 0.02;
            PbodyHydro(count)->ZIa = 0.02;
#if (UseSFModelSpawn) 
            PbodyHydro(count)->SpawnMass = mgas/((double)MaxSpawnTimes);
#endif
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
            PbodyHydro(count)->G0 = 10;
            PbodyHydro(count)->fH2 = 0.4;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

            PbodyHydro(count)->DQheat = 0.e0;
            PbodyHydro(count)->Use = ON;
            PbodyHydro(count)->Kernel = eps;
            PbodyHydro(count)->U = Uinit;
#endif
            count ++;
        }
    }
#else
    while(count < Pall.Nhydro){
    //while(count < Pall.Ntotal){
        double Pos[3];
        Pos[0] = Rmax*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        Pos[1] = Rmax*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        double r = sqrt(SQ(Pos[0])+SQ(Pos[1]));
        if((r>1*PC_CGS/Pall.UnitLength) && (r<Rmax)){

            Pbody[count]->Active = ON;
            Pbody[count]->Use = ON;
            Pbody[count]->Type = TypeHydro;
            //Pbody[count]->Type = TypeDM;
            Pbody[count]->GlobalID = count+Pall.Ntotal*MyID;

            Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
            Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
            Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Zsig*Gaussian();

            //double Vc = CircularVelocityFromExternalPotentials(r);
            double Vc = InitializeCircularVelocityFromAGNExternalPotentials(count,r);
            Pbody[count]->Vel[0] = -Vc*Pos[1]/r;
            Pbody[count]->Vel[1] =  Vc*Pos[0]/r;
            Pbody[count]->Vel[2] = 0.e0;

            Pbody[count]->Vel[0] += 0.3*Vc*Gaussian();
            Pbody[count]->Vel[1] += 0.3*Vc*Gaussian();
            Pbody[count]->Vel[2] += 0.3*Vc*Gaussian();

            // Pbody[count]->Vel[0] += cs*Gaussian();
            // Pbody[count]->Vel[1] += cs*Gaussian();
            // Pbody[count]->Vel[2] += cs*Gaussian();

            Pbody[count]->Eps = eps;
            Pbody[count]->Mass = mgas;
#if 1
            PbodyHydro(count)->Mass = Pbody[count]->Mass = mgas;
            PbodyHydro(count)->Z = 0.02;
#if (UseSFModelSpawn) 
            PbodyHydro(count)->SpawnMass = mgas/((double)MaxSpawnTimes);
#endif
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
            PbodyHydro(count)->G0 = 10;
            PbodyHydro(count)->fH2 = 0.4;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

            PbodyHydro(count)->Use = ON;
            PbodyHydro(count)->Kernel = eps;
            PbodyHydro(count)->U = Uinit;
            // if(gsl_rng_uniform(RandomGenerator) < 0.01)
                // PbodyHydro(count)->U = 1.e+3*Uinit;
#endif
            count ++;
        } else {
        }
    }
#endif

    if(MPIGetMyID() == MPI_ROOT_RANK){
        int index = Pall.Ntotal-1;
        Pbody[index]->Active = ON;
        Pbody[index]->Use = ON;
        Pbody[index]->Type = TypeSink;
        Pbody[index]->GlobalID = Pall.Ntotal_t-1;

        Pbody[index]->Pos[0] = Pbody[index]->Pos[1] = Pbody[index]->Pos[2] = 0.0;
        Pbody[index]->Vel[0] = Pbody[index]->Vel[1] = Pbody[index]->Vel[2] = 0.e0;

        Pbody[index]->Eps = eps;
        Pbody[index]->Mass = MassBH_AGN_CGS/Pall.UnitMass;
        //Pbody[index]->Mass = MBH_AGN/Pall.UnitMass;

        PbodySink(index)->Use = ON;
        PbodySink(index)->ParentGlobalID = Pbody[index]->GlobalID;
        PbodySink(index)->FormationTime = 0.e0;
        PbodySink(index)->PosP[0] = PbodySink(index)->PosP[1] = PbodySink(index)->PosP[2] = 0.e0;
        PbodySink(index)->VelP[0] = PbodySink(index)->VelP[1] = PbodySink(index)->VelP[2] = 0.e0;

        PbodySink(index)->AccretionRadius = SINKHYDRO_ACCRETION_RADIUS/Pall.UnitLength;
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
#ifndef __AGN_TORUS_TEST__ //FB test
        PhydroBody(i)->Vel[2] = cs*Gaussian(); //FB test
#endif
    }
    
    for(int i=0;i<Pall.Ntotal;i++){
        Pbody[i]->Active = ON;
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;
    }
    AGNForceFromExternalPotentials();

#define VSigma (0.1)
    double Acc_gas = Pall.GravConst*(MassGas_AGN_CGS/Pall.UnitMass)/SQ(Rmax);
    //double Acc_gas = Pall.GravConst*(MGas_AGN/Pall.UnitMass)/SQ(Rmax);
    for(int i=0;i<Pall.Ntotal;i++){
        double R = sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1]));
        // Selfgravity of gaseous disk.
        if(R > 0.e0){
#ifdef USE_HERNQUIST_BULGE_FOR_AGN
            double vc_Gas = VcGas_AGN(R);
            double vc_BH = VcBH_AGN(R,eps);
            double vc_Bulge = VcBulge_AGN(R);
            double vc = sqrt(SQ(vc_Gas)+SQ(vc_BH)+SQ(vc_Bulge));
            //double vc = sqrt(SQ(vc_Gas));
            //double vc = sqrt(SQ(vc_BH));
            // double vc = sqrt(SQ(vc_Bulge));
            if(R>1.0*PC_CGS/Pall.UnitLength){
                Pbody[i]->Vel[0] = -(vc/R)*Pbody[i]->Pos[1];
                Pbody[i]->Vel[1] = +(vc/R)*Pbody[i]->Pos[0];
            }  else {
                Pbody[i]->Vel[0] = 1.e-10*GetUnitVel();
                Pbody[i]->Vel[1] = 1.e-10*GetUnitVel();
            }

#if 0
            // Gas disk contribution.
            Pbody[i]->Acc[0] -= Acc_gas*Pbody[i]->Pos[0]/R; // contribution from gas
            Pbody[i]->Acc[1] -= Acc_gas*Pbody[i]->Pos[1]/R; // contribution from gas
            Pbody[i]->Acc[2] = 0.0;                         // contribution from gas
#endif

#if 0
            // BH contribution. 
            double r = sqrt(R*R+eps*eps);
            double Acc_BH = Pall.GravConst*(MBH_AGN/Pall.UnitMass)/CUBE(r);
            Pbody[i]->Acc[0] -= Acc_BH*Pbody[i]->Pos[0];   // contribution from BH
            Pbody[i]->Acc[1] -= Acc_BH*Pbody[i]->Pos[1];   // contribution from BH
            Pbody[i]->Acc[2] -= Acc_BH*Pbody[i]->Pos[2];   // contribution from BH
#endif

#if 0
            // Bulge contribution.
            double rhern = (R+RaHernquistUnit);
            double Acc_Bulge = GMBulgeUnit/CUBE(rhern);
            Pbody[i]->Acc[0] -= Acc_Bulge*Pbody[i]->Pos[0];   // contribution from Bulge
            Pbody[i]->Acc[1] -= Acc_Bulge*Pbody[i]->Pos[1];   // contribution from Bulge
            Pbody[i]->Acc[2] -= Acc_Bulge*Pbody[i]->Pos[2];   // contribution from Bulge
#endif
#else
            // Pbody[i]->Acc[0] -= Acc_gas*Pbody[i]->Pos[0]/R;
            //Pbody[i]->Acc[1] -= Acc_gas*Pbody[i]->Pos[1]/R;
            Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.0;

            double r = sqrt(R*R+eps*eps);
            // double Acc_BH = Pall.GravConst*(MBH_AGN/Pall.UnitMass)/CUBE(r);
            double ra1 = sqrt(R*R+SQ(Ra1_AGN/Pall.UnitLength));
            double ra2 = sqrt(R*R+SQ(Ra2_AGN/Pall.UnitLength));
            //double Acc_BH = Pall.GravConst*(MBH_AGN/Pall.UnitMass)/CUBE(r);
            // double Acc_BH = Pall.GravConst*(MBH_AGN/Pall.UnitMass)/CUBE(r)+
                            // FVcUnit2*(1.0/CUBE(ra1)+1.0/CUBE(ra2));
            double Acc_BH = FVcUnit2*(1.0/CUBE(ra1)+1.0/CUBE(ra2));
            Pbody[i]->Acc[0] -= Acc_BH*Pbody[i]->Pos[0];
            Pbody[i]->Acc[1] -= Acc_BH*Pbody[i]->Pos[1];
            Pbody[i]->Acc[2] -= Acc_BH*Pbody[i]->Pos[2];

            double Acc = sqrt(SQ(Pbody[i]->Acc[0])+SQ(Pbody[i]->Acc[1]));
            Pbody[i]->Vel[0] = -sqrt(Acc/R)*Pbody[i]->Pos[1];
            Pbody[i]->Vel[1] = +sqrt(Acc/R)*Pbody[i]->Pos[0];
#endif

#ifndef __AGN_TORUS_TEST__ //FB test
            if(Pbody[i]->Type == TypeHydro)
                Pbody[i]->Vel[2] = gsl_ran_gaussian(RandomGenerator,sqrt(5.0/3.0*(2.0/3.0)*PbodyHydro(i)->U));
#endif

            double factor_vc = sqrt(SQ(Pbody[i]->Vel[0])+SQ(Pbody[i]->Vel[1]));

#ifndef __AGN_TORUS_TEST__ //FB test
            Pbody[i]->Vel[0] *= (1.0+gsl_ran_gaussian(RandomGenerator,VSigma));
            Pbody[i]->Vel[1] *= (1.0+gsl_ran_gaussian(RandomGenerator,VSigma));
#endif

            double factor_vcvs = sqrt(SQ(Pbody[i]->Vel[0])+SQ(Pbody[i]->Vel[1]));
            Pbody[i]->Vel[0] *= factor_vc/factor_vcvs;
            Pbody[i]->Vel[1] *= factor_vc/factor_vcvs;
        }
        //fprintf(stderr,"%g %g %g %g\n",Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Acc[0],Pbody[i]->Acc[1]);
    }
    //fflush(stderr);
    //exit(1);
    //
    
#if 0
    {

    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"./init.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]);
    }
    fclose(fp_init);
    fflush(NULL);
    MPI_Finalize();
    exit(1); 
    }
#endif
 

    if(DistortionFlag == 0){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"DistotionFlag = %d\n",DistortionFlag);
            fprintf(stderr,"Nothing is done.\n");
        }
    } else if(DistortionFlag == 1){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"DistotionFlag = %d\n",DistortionFlag);
            fprintf(stderr,"Stretch particle dist. five times, uniformly.\n");
        }
        for(int i=0;i<Pall.Nhydro;i++){
            double fr = 5;
            double fv = 1.0/fr;
            PhydroBody(i)->Pos[0] *= fr;
            PhydroBody(i)->Pos[1] *= fr;

            PhydroBody(i)->Vel[0] *= fv;
            PhydroBody(i)->Vel[1] *= fv;
        }
    } else if(DistortionFlag == 2){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"DistotionFlag = %d\n",DistortionFlag);
            fprintf(stderr,"Stretch particle dist. following this mapping function: f(r) = (1+4*r/R).\n");
        }
        for(int i=0;i<Pall.Nhydro;i++){
            double r = sqrt(SQ(PhydroBody(i)->Pos[0])+SQ(PhydroBody(i)->Pos[1]));
            double fr = (1.0+4*r/Rmax);
            double fv = 1.0/fr;
            PhydroBody(i)->Pos[0] *= fr;
            PhydroBody(i)->Pos[1] *= fr;

            PhydroBody(i)->Vel[0] *= fv;
            PhydroBody(i)->Vel[1] *= fv;
        }
    } else if(DistortionFlag == 3){
        if(MPIGetMyID() == MPI_ROOT_RANK){
            fprintf(stderr,"DistotionFlag = %d\n",DistortionFlag);
            fprintf(stderr,"Stretch and distort particle distribution.\n");
        }
        // AM map function.
        double StretchFact = 5;
        double InnerAM = 0.2;
        double OuterAM = 1.0;
        double DiskInEdge = (1*PC_CGS/Pall.UnitLength)*StretchFact;
        double DiskOutEdge = (Redge_AGN_CGS/Pall.UnitLength)*StretchFact;

        double CoefA = (OuterAM-InnerAM)/DiskOutEdge;
        double CoefB = InnerAM;

        fprintf(stderr,"F,IAM,OAM,DI,DE,CA,CB = %g %g %g %g %g %g %g\n",
                StretchFact,InnerAM,OuterAM,DiskInEdge,DiskOutEdge,CoefA,CoefB);

#if 0
        FILE *fp;
        FileOpen(fp,"l.data","w");
        for(int i=0;i<Pall.Nhydro;i++){
            double r = sqrt(SQ(PhydroBody(i)->Pos[0])+SQ(PhydroBody(i)->Pos[1]));
            double lz = 
                (PhydroBody(i)->Pos[0])*(PhydroBody(i)->Vel[1])-
                (PhydroBody(i)->Pos[1])*(PhydroBody(i)->Vel[0]);
            fprintf(fp,"%g %g\n",r,lz);
        }
        fclose(fp);
        fflush(fp);
        exit(1);
#endif

        // FILE *fp;
        // FileOpen(fp,"l.data","w");
        for(int i=0;i<Pall.Nhydro;i++){
            Pbody[i]->Pos[0] *= StretchFact;
            Pbody[i]->Pos[1] *= StretchFact;
            double R = sqrt(SQ(Pbody[i]->Pos[0])+SQ(Pbody[i]->Pos[1]));
            double vc_Gas = VcGas_AGN(R);
            double vc_BH = VcBH_AGN(R,eps);
            double vc_Bulge = VcBulge_AGN(R);
            double vc = sqrt(SQ(vc_Gas)+SQ(vc_BH)+SQ(vc_Bulge));
            if(R>1.0*PC_CGS/Pall.UnitLength){
                Pbody[i]->Vel[0] = -(vc/R)*Pbody[i]->Pos[1];
                Pbody[i]->Vel[1] = +(vc/R)*Pbody[i]->Pos[0];
            }  else {
                Pbody[i]->Vel[0] = 1.e-10*GetUnitVel();
                Pbody[i]->Vel[1] = 1.e-10*GetUnitVel();
            }

            // double r = sqrt(SQ(PhydroBody(i)->Pos[0])+SQ(PhydroBody(i)->Pos[1]));
            // double lz = 
                // (PhydroBody(i)->Pos[0])*(PhydroBody(i)->Vel[1])-
                // (PhydroBody(i)->Pos[1])*(PhydroBody(i)->Vel[0]);


            double ReduceFact = CoefA*R+CoefB;
            Pbody[i]->Vel[0] *= ReduceFact;
            Pbody[i]->Vel[1] *= ReduceFact;
            // double lz2 = 
                // (PhydroBody(i)->Pos[0])*(PhydroBody(i)->Vel[1])-
                // (PhydroBody(i)->Pos[1])*(PhydroBody(i)->Vel[0]);
            // fprintf(fp,"%g %g %g\n",r,lz,lz2);
        }
        // fclose(fp);
        // fflush(fp);
        // exit(1);

    } else {
        if(MPIGetMyID() == MPI_ROOT_RANK)
            fprintf(stderr,"Need correct setting of DistotionFlag\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(1);
    }

#if 0
#if 1 // particle stretch type 1: f(r) = (1+4*r/R) map.
    for(int i=0;i<Pall.Nhydro;i++){
        double r = sqrt(SQ(PhydroBody(i)->Pos[0])+SQ(PhydroBody(i)->Pos[1]));
        double fr = (1.0+4*r/Rmax);
        double fv = 1.0/fr;
        PhydroBody(i)->Pos[0] *= fr;
        PhydroBody(i)->Pos[1] *= fr;

        PhydroBody(i)->Vel[0] *= fv;
        PhydroBody(i)->Vel[1] *= fv;
    }
#else // particle stretch type 2: constant separation. r x v = constant.
    double r_new = 5.0;
    double v_new = 1.0/r_new;
    for(int i=0;i<Pall.Nhydro;i++){
        PhydroBody(i)->Pos[0] *= r_new;
        PhydroBody(i)->Pos[1] *= r_new;

        PhydroBody(i)->Vel[0] *= v_new;
        PhydroBody(i)->Vel[1] *= v_new;
    }
#endif
#endif

#undef __AGN_TORUS_TEST__

    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"./init.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]);
    }
    fclose(fp_init);

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    //
    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.TEnd = 100*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
                Pall.TEnd,Pall.TEnd*Pall.UnitTime);
    }


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.0001*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/AGN.ASCII");
    strcpy(Pall.BaseFileName,"./data/AGN");
    strcpy(Pall.RestartFileName,"./data/AGN.dump");

    return;
}
#endif //TASK_AGN


#ifdef TASK_GALACTIC_CENTER
extern int GCCloudType;
extern double GCCloudPower;
extern int GCCloudN;
extern double GCCloudR;
extern double GCCloudM;
extern double GCMBH;
extern double GCEPSBH;
extern double GCEATHALOGAS;
extern double GCR_peri;
extern double GCDay_peri;
extern double GCEccentricity;
extern double GCTorbit;
extern double GCDay_start;
extern double GCDay_end;
extern int GCNorbit;
extern char GCOutputDir[MaxCharactersInLine];
// halo gas parameter
extern int GCHaloType; // Type 0 simple hydrostatic
extern int GCHaloN;
extern int GCHaloRotation; // if 1 halo has angular momentum.
                           // This flag works only when GCHaloType == 0
extern double GCHaloRotationFactor; // if GCHaloRotation == 1, halo rotates with ThisParam X keplar.
extern int GCHaloSpinVectorChange; // if 1, change the spin vector's direction.
extern double GCHaloSpinTheta;
extern double GCHaloSpinPhi;
extern double GCHaloEtaHot;
extern double GCHaloTruncR;
// halo gas parameter

static void GCSetInitParticleData(const double Pos[], const int Index, const int GlobalID, const double eps, const double mgas, const double Uinit, const int Tag){

    for(int k=0;k<3;k++)
        Pbody[Index]->PosP[k] = Pbody[Index]->Pos[k] = Pos[k];

    Pbody[Index]->Active = ON;
    Pbody[Index]->Use = ON;
    Pbody[Index]->Type = TypeHydro;
    Pbody[Index]->GlobalID = GlobalID;

    Pbody[Index]->Eps = eps;
    Pbody[Index]->Mass = mgas;

    PbodyHydro(Index)->Mass = Pbody[Index]->Mass = mgas;
    PbodyHydro(Index)->Z = 0.02;
    // PbodyHydro(Index)->ZII = 0.02;
    // PbodyHydro(Index)->ZIa = 0.02;
    PbodyHydro(Index)->ZII = 
    PbodyHydro(Index)->ZIa = 0.0;
#if (UseSFModelSpawn) 
    PbodyHydro(Index)->SpawnMass = mgas/((double)MaxSpawnTimes);
#endif
#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
    PbodyHydro(Index)->G0 = 1.7e7;
    PbodyHydro(Index)->fH2 = 0.1;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

    PbodyHydro(Index)->DQheat = 0.e0;
    PbodyHydro(Index)->Use = ON;
    PbodyHydro(Index)->Kernel = eps;
    PbodyHydro(Index)->U = Uinit;

#ifdef USE_PARTICLE_TAG
    PbodyHydro(Index)->Tag = Tag;
#endif

    return ;
}

static double GCReturnCloudMass(double Radius){

    if((GCCloudPower > -3)||(GCCloudPower < -3)){
        return (GCCloudM/pow(GCCloudR,3.0+GCCloudPower))*pow(Radius,3.0+GCCloudPower);
    } else {
        return (GCCloudM/log10(GCCloudR))*log10(Radius);
    }
}

static double GCBSearchNextRadiusCloud(const double Mass, const double PrevRadius){

    double RadiusNext;
    //double RadiusLeft = 0.e0;
    double RadiusLeft = fmax(0.e0,PrevRadius);
    double RadiusRight = 2.0*GCCloudR;

    int Iteration = 0;
    do{
        RadiusNext = 0.5*(RadiusLeft + RadiusRight);
        double CurrentMass = GCReturnCloudMass(RadiusNext);
        if(CurrentMass == Mass) return RadiusNext;
        // gprintlmpi(CurrentMass);
        // gprintlmpi(RadiusNext);
        if(CurrentMass >  Mass){
            RadiusRight = RadiusNext;
        }
        if(CurrentMass <  Mass){
            RadiusLeft = RadiusNext;
        }
        // gprintlmpi(fabs(RadiusLeft-RadiusRight));
        Iteration++;
        if(Iteration >100) exit(1);
    } while (fabs(RadiusLeft-RadiusRight)>1.e-6*RadiusRight);

    return RadiusNext;
}


static double GCHaloRho0 = 1.7e-21;
static double GCHaloR0 = 1.e+16;
static double GCReturnHaloMass(const double R){ // unit in gram.
    if(GCHaloType == 0){
        return GCHaloEtaHot*2.0*M_PI*GCHaloRho0*GCHaloR0*SQ(R*Pall.UnitLength)/Pall.UnitMass;
    } else {
        return -1;
    }
}
static double GCReturnHaloTemperature(const double R){ // unit in Kelvin.
    if(GCHaloType == 0){
        return 2.1e+8*(GCHaloR0/(R*Pall.UnitLength));
        //return 3.5e+8*(GCHaloR0/(R*Pall.UnitLength));
    } else { 
        return -1;
    }
}

static double GCBSearchNextRadiusHalo(const double Mass, const double PrevRadius){

    double RadiusNext;
    //double RadiusLeft = 0.e0;
    double RadiusLeft = fmax(0.e0,PrevRadius);
    double RadiusRight = 2.0*GCHaloTruncR;

    int Iteration = 0;
    do{
        RadiusNext = 0.5*(RadiusLeft + RadiusRight);
        double CurrentMass = GCReturnHaloMass(RadiusNext);
        if(CurrentMass == Mass) return RadiusNext;
        if(CurrentMass >  Mass){
            RadiusRight = RadiusNext;
        }
        if(CurrentMass <  Mass){
            RadiusLeft = RadiusNext;
        }
        Iteration++;
        if(Iteration >100) exit(1);
    } while (fabs(RadiusLeft-RadiusRight)>1.e-6*RadiusRight);

    return RadiusNext;
}

static void GCWriteHaloProfile(const int Number, const int NumberCloud, const int NumberHalo){ // Write halo profile

    int Nbin = 100;
    double dLogR = (log10(GCHaloTruncR) - log10(GCHaloTruncR/1000.0))/(double)Nbin;
    gprintlmpi(log10(GCHaloTruncR));
    gprintlmpi(log10(GCHaloTruncR/1000.0));
    gprintlmpi(dLogR);

    double MassBin[Nbin];
    double DensityBin[Nbin];
    double TemperatureBin[Nbin];
    for(int i=0;i<Nbin;i++){
        MassBin[i] = 0;
        TemperatureBin[i] = 0;
    }
   
    for(int i=NumberCloud;i<Number;i++){
        double LogR = log10(NORM(PhydroBody(i)->Pos))-log10(GCHaloTruncR/1000.0);
        int IntLogR = MAX(0,LogR/dLogR);
        // fprintf(stderr,"%g %g %d\n",LogR,dLogR,IntLogR);
        if(IntLogR >= Nbin) continue;
        assert(IntLogR >= 0);
        MassBin[IntLogR] += Phydro[i]->Mass;
        TemperatureBin[IntLogR] += Phydro[i]->Mass*Phydro[i]->U;
    }
    double GlobalMassBin[Nbin];
    double GlobalTemperatureBin[Nbin];
    MPI_Allreduce(MassBin,GlobalMassBin,Nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(TemperatureBin,GlobalTemperatureBin,Nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for(int i=0;i<Nbin;i++){
        double Rin  = pow(10.0,dLogR*i);
        double Rout = pow(10.0,dLogR*(i+1));
        double dV = (4*M_PI/3.0)*(CUBE(Rout)-CUBE(Rin));
        if(GlobalMassBin[i] > 0.e0){
            GlobalTemperatureBin[i] /= GlobalMassBin[i];
            GlobalTemperatureBin[i] *= Pall.ConvertUtoT;

            GlobalMassBin[i]/=dV;
        }
    }


    if(MPIGetMyID() == MPI_ROOT_RANK){
        FILE *fp;
        char Fname[MaxCharactersInLine];
        Snprintf(Fname,"%s/InitialHaloProfile.dat",GCOutputDir);
        FileOpen(fp,Fname,"w");
        fprintf(fp,"#This_routine_would_be_wrong.");
        for(int i=0;i<Nbin;i++){
            double Rmean  = pow(10.0,dLogR*(i+0.5));
            fprintf(fp,"%g %g %g\n",Rmean,GlobalMassBin[i],GlobalTemperatureBin[i]);
        }
        fclose(fp);
    }

    return ;
}

double GCInitPosForRemove[3];
static bool GCParticleRemoveCondition(const int index){

    if(Phydro[index]->Tag == 1){
        if(DISTANCE(PhydroBody(index)->Pos,GCInitPosForRemove) < GCCloudR){
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

void InitGCCloud(double InitPos[], double InitVel[]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitLength = AU_CGS;
    Pall.UnitTime = YEAR_CGS;
    Pall.UnitMass = MEARTH_CGS;

    // Convert Parameters into simulation units.
    GCCloudR /= Pall.UnitLength;
    GCCloudM /= Pall.UnitMass;
    GCMBH    /= Pall.UnitMass;
    GCEPSBH  /= Pall.UnitLength;
    GCR_peri /= Pall.UnitLength;
    GCHaloTruncR /= Pall.UnitLength;

    fprintf(stderr,"%g %g %g\n",
            InitPos[0],InitPos[1],InitPos[2]);
    for(int k=0;k<3;k++){
        InitPos[k] /= Pall.UnitLength;
        InitVel[k] *= Pall.UnitTime/Pall.UnitLength;
    }
    fprintf(stderr,"%g %g %g\n",
            InitPos[0],InitPos[1],InitPos[2]);

    // Convert Parameters into simulation units.
    // Convert orbit infomation.
    if(MyID == MPI_ROOT_RANK){
        FILE *fp_input,*fp_output;
        char fname_input[MaxCharactersInLine];
        char fname_output[MaxCharactersInLine];
        Snprintf(fname_input,"%s/Orbit.dat",GCOutputDir);
        FileOpen(fp_input,fname_input,"r");
        Snprintf(fname_output,"%s/OrbitUnit.dat",GCOutputDir);
        FileOpen(fp_output,fname_output,"w");
        double PosOrbit[3];
        while(fscanf(fp_input,"%le %le %le",PosOrbit, PosOrbit+1, PosOrbit+2)!=EOF){
            PosOrbit[0] /= Pall.UnitLength;
            PosOrbit[1] /= Pall.UnitLength;
            PosOrbit[2] /= Pall.UnitLength;
            fprintf(fp_output,"%g %g %g\n",PosOrbit[0],PosOrbit[1],PosOrbit[2]);
        }
        fclose(fp_input);
        fclose(fp_output);


        Snprintf(fname_input,"%s/Integral.log",GCOutputDir);
        FileOpen(fp_input,fname_input,"r");
        Snprintf(fname_output,"%s/IntegralUnit.log",GCOutputDir);
        FileOpen(fp_output,fname_output,"w");
        double VelOrbit[3];
        double AccOrbit[3];
        while(fscanf(fp_input,"%le %le %le %le %le %le %le %le %le",
                    PosOrbit, PosOrbit+1, PosOrbit+2,
                    VelOrbit, VelOrbit+1, VelOrbit+2,
                    AccOrbit, AccOrbit+1, AccOrbit+2)!=EOF){
            PosOrbit[0] /= Pall.UnitLength;
            PosOrbit[1] /= Pall.UnitLength;
            PosOrbit[2] /= Pall.UnitLength;

            VelOrbit[0] *= Pall.UnitTime/Pall.UnitLength;
            VelOrbit[1] *= Pall.UnitTime/Pall.UnitLength;
            VelOrbit[2] *= Pall.UnitTime/Pall.UnitLength;

            AccOrbit[0] *= SQ(Pall.UnitTime)/Pall.UnitLength;
            AccOrbit[1] *= SQ(Pall.UnitTime)/Pall.UnitLength;
            AccOrbit[2] *= SQ(Pall.UnitTime)/Pall.UnitLength;
            fprintf(fp_output,"%g %g %g %g %g %g %g %g %g\n",
                    PosOrbit[0],PosOrbit[1],PosOrbit[2],
                    VelOrbit[0],VelOrbit[1],VelOrbit[2],
                    AccOrbit[0],AccOrbit[1],AccOrbit[2]);
        }
        fclose(fp_input);
        fclose(fp_output);
    }
    // Convert orbit infomation.

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    int NumberCloud = GCCloudN;
    int NumberHalo = GCHaloN;
    int Number = NumberCloud+NumberHalo;

    Pall.Nhydro = NumberCloud/NProcs+NumberHalo/NProcs;
    Pall.Nhydro_t = Pall.Nhydro*NProcs;
    if(MyID == MPI_ROOT_RANK){
        Pall.Nsink = 1;
    } else {
        Pall.Nsink = 0;
    }
    Pall.Nsink_t = 1;
    Pall.Ntotal = Pall.Nhydro + Pall.Nsink;
    Pall.Ntotal_t = Pall.Nhydro_t + Pall.Nsink_t;

    NumberCloud /= NProcs;
    NumberHalo /= NProcs;
    Number = NumberCloud+NumberHalo;
    int NumberCloud_t = NumberCloud*NProcs;
    int NumberHalo_t = NumberHalo*NProcs;
    
    InitializeAllActiveParticleNumbers();

    int AllocationSize = Pall.Ntotal; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    if(MyID == MPI_ROOT_RANK){
        GenerateStructPsink(1);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    if(MyID == MPI_ROOT_RANK){
        Pbody[AllocationSize-1]->Baryon = (void *)(Psink[0]);
        Psink[0]->Body = Pbody[AllocationSize-1];
    }
    int SinkID = Pall.Ntotal-1;
    
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,6.0,1.0,0.1);
#else
    SetConstantViscosityParameters(6.0);
#endif // USE_VARIABLE_ALPHA //}

    double mgas = GCCloudM/MAX((double)NumberCloud_t,1);
    fprintf(stderr,"%g [MEarth], NCloud %d\n",mgas,NumberCloud_t);
    
    double eps = cbrt(0.01*(4.0*M_PI/3.0)*CUBE(GCCloudR)/MAX((double)NumberCloud_t,1));
    fprintf(stderr,"Eps = %g AU\n",eps*Pall.UnitLength/AU_CGS);
    double Uinit = 10000*Pall.ConvertTtoU;
    double cs = sqrt(Pall.GGm1*Uinit);
    fprintf(stderr,"Sound speed = %g km/s\n",1.e-5*(cs*Pall.UnitLength/Pall.UnitTime));

#if 0
    if(MyID == MPI_ROOT_RANK){
        if(NumberCloud_t != 0)
            GCSetInitParticleData((double[]){0.0,0.0,0.0},0,0,eps,mgas,Uinit,0);
    }

    int counter;
    if(MyID == MPI_ROOT_RANK){
        counter = 1;
    } else {
        counter = 0;
    }
    // FILE *ff;
    // FileOpen(ff,"./aa","w")
    double CurrentMass = mgas;
    for(int i=1;i<NumberCloud_t;i++){
        double Rand[2] = {gsl_rng_uniform(RandomGenerator),
                          gsl_rng_uniform(RandomGenerator)};
        double PP;
        if(i%NProcs ==  MyID){
            // gprintlmpi(CurrentMass);
            double Radius = GCBSearchNextRadiusCloud(CurrentMass);
            // PP = Radius;
            // gprintlmpi(Radius);
            double phi = 2*M_PI*Rand[0];
            double cos_theta = (2.0*Rand[1]-1.0);
            double sin_theta = sqrt(1.0-SQ(cos_theta));
            double Pos[3] = {Radius*sin_theta*cos(phi),
                             Radius*sin_theta*sin(phi),
                             Radius*cos_theta};
            // fprintf(stderr,"R %g %g, Pos %g %g %g\n",Rand[0],Rand[1],Pos[0],Pos[1],Pos[2]);
            GCSetInitParticleData(Pos,counter,i,eps,mgas,Uinit,0);
            // gprintlmpi(Pos[2]);
            counter ++;
        }
        CurrentMass = (i+1)*mgas;
        // fprintf(ff,"%g %g\n",PP,CurrentMass);
    }
    // fclose(ff);
#else
    int counter = 0;
    double PrevRadius = 0.e0;
    for(int i=0;i<NumberCloud_t;i++){
        double CurrentMass = (i+1)*mgas;
        double Rand[2] = {gsl_rng_uniform(RandomGenerator),
                          gsl_rng_uniform(RandomGenerator)};
        double PP;
        if(i%NProcs ==  MyID){
            double Radius = GCBSearchNextRadiusCloud(CurrentMass,PrevRadius);
            double phi = 2*M_PI*Rand[0];
            double cos_theta = (2.0*Rand[1]-1.0);
            double sin_theta = sqrt(1.0-SQ(cos_theta));
            double Pos[3] = {Radius*sin_theta*cos(phi),
                             Radius*sin_theta*sin(phi),
                             Radius*cos_theta};
            GCSetInitParticleData(Pos,counter,i,eps,mgas,Uinit,0);
            counter ++;
            PrevRadius = Radius;
        }
    }
#endif

    // Halo loop
    double mhalo = GCReturnHaloMass(GCHaloTruncR)/MAX((double)NumberHalo_t,1);
    double epshalo;
    if(GCCloudN > 0){
        epshalo = cbrt(mhalo/mgas)*eps;
    } else {
        double _mgas = GCCloudM/(double)1000;
        double _eps = cbrt(0.01*(4.0*M_PI/3.0)*CUBE(GCCloudR)/(double)1000);
        epshalo = cbrt(mhalo/_mgas)*_eps;
    }
    fprintf(stderr,"Halo Mass = %g [EM], mhalo = %g [EM], eps_halo = %g [AU]\n",
            GCReturnHaloMass(GCHaloTruncR),mhalo,epshalo);


#if 0
    double RandHalo[2] = {gsl_rng_uniform(RandomGenerator),
                          gsl_rng_uniform(RandomGenerator)};
    if(MyID == MPI_ROOT_RANK){
        if(NumberHalo_t != 0){
            double Radius = GCBSearchNextRadiusHalo(mhalo);
            double phi = 2*M_PI*RandHalo[0];
            double cos_theta = (2.0*RandHalo[1]-1.0);
            double sin_theta = sqrt(1.0-SQ(cos_theta));
            double Pos[3] = {Radius*sin_theta*cos(phi),
                             Radius*sin_theta*sin(phi),
                             Radius*cos_theta};
            double Uhalo = GCReturnHaloTemperature(Radius)*Pall.ConvertTtoU;
            GCSetInitParticleData(Pos,counter,NumberCloud_t,epshalo,mhalo,Uhalo,1);
            counter ++;
        }
    }
    CurrentMass = 2*mhalo;
    for(int i=1;i<NumberHalo_t;i++){ //loop for halo
        double RandHalo[2] = {gsl_rng_uniform(RandomGenerator),
                              gsl_rng_uniform(RandomGenerator)};
        if(i%NProcs ==  MyID){
            double Radius = GCBSearchNextRadiusHalo(CurrentMass);
            double phi = 2*M_PI*RandHalo[0];
            double cos_theta = (2.0*RandHalo[1]-1.0);
            double sin_theta = sqrt(1.0-SQ(cos_theta));
            double Pos[3] = {Radius*sin_theta*cos(phi),
                             Radius*sin_theta*sin(phi),
                             Radius*cos_theta};
            /// 
            double Uhalo = GCReturnHaloTemperature(Radius)*Pall.ConvertTtoU;
            GCSetInitParticleData(Pos,counter,NumberCloud_t+i,epshalo,mhalo,Uhalo,1);
            counter ++;
        }
        CurrentMass = (i+2)*mhalo;
    }
#else
    PrevRadius = 0.e0;
    for(int i=0;i<NumberHalo_t;i++){ //loop for halo
        double CurrentMass = (i+1)*mhalo;
        double RandHalo[2] = {gsl_rng_uniform(RandomGenerator),
                              gsl_rng_uniform(RandomGenerator)};
        double Radius = GCBSearchNextRadiusHalo(CurrentMass,PrevRadius);
        double phi = 2*M_PI*RandHalo[0];
        double cos_theta = (2.0*RandHalo[1]-1.0);
        double sin_theta = sqrt(1.0-SQ(cos_theta));
        double Pos[3] = {Radius*sin_theta*cos(phi),
                         Radius*sin_theta*sin(phi),
                         Radius*cos_theta};
        if(i%NProcs ==  MyID){
            // if((NumberCloud_t>0)&&(DISTANCE(Pos,InitPos) < GCCloudR)){
                // Pall.Nhydro --;
                // Pall.Ntotal --;
            // } else {
                double Uhalo = GCReturnHaloTemperature(Radius)*Pall.ConvertTtoU;
                GCSetInitParticleData(Pos,counter,NumberCloud_t+i,epshalo,mhalo,Uhalo,1);
                counter ++;
            // }
            PrevRadius = Radius;
        }
    }
#endif

    //double eps_BH = eps;
    if(MyID == MPI_ROOT_RANK){
        int index = SinkID;
        Pbody[index]->Active = ON;
        Pbody[index]->Use = ON;
        Pbody[index]->Type = TypeSink;
        Pbody[index]->GlobalID = Pall.Ntotal_t-1;

        Pbody[index]->Pos[0] = Pbody[index]->Pos[1] = Pbody[index]->Pos[2] = 0.0;
        Pbody[index]->Vel[0] = Pbody[index]->Vel[1] = Pbody[index]->Vel[2] = 0.e0;

        Pbody[index]->Eps = GCEPSBH;
        Pbody[index]->Mass = GCMBH;

        PbodySink(index)->Use = ON;
        PbodySink(index)->ParentGlobalID = Pbody[index]->GlobalID;
        PbodySink(index)->FormationTime = 0.e0;
        PbodySink(index)->PosP[0] = PbodySink(index)->PosP[1] = PbodySink(index)->PosP[2] = 0.e0;
        PbodySink(index)->VelP[0] = PbodySink(index)->VelP[1] = PbodySink(index)->VelP[2] = 0.e0;

        PbodySink(index)->AccretionRadius = SINKHYDRO_ACCRETION_RADIUS/Pall.UnitLength;
    }

    ////////////////////////////////////////////////////////////////////////////


    // Shift Positions and Velocities.
    for(int i=0;i<Pall.Nhydro;i++){
        for(int k=0;k<3;k++){
            if(i < NumberCloud){
                PhydroBody(i)->Pos[k] += InitPos[k];
                PhydroBody(i)->PosP[k] += InitPos[k];
                PhydroBody(i)->Vel[k] += InitVel[k];
            }
            Phydro[i]->PosP[k] = PhydroBody(i)->Pos[k];
        }
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

    if(GCHaloRotation != 0){
        for(int i=0;i<Pall.Nhydro;i++){
            if(Phydro[i]->Tag == 1){
                double r = sqrt(NORM2(PhydroBody(i)->Pos)+SQ(GCEPSBH));
                // double v_tan = GCHaloRotationFactor*sqrt(Pall.GravConst*GCMBH/r);
                double v_tan = GCHaloRotationFactor*sqrt(Pall.GravConst*GCMBH/r
                        *(sqrt(SQ(PhydroBody(i)->Pos[0])+SQ(PhydroBody(i)->Pos[1]))/NORM(PhydroBody(i)->Pos)));
                double theta = atan(PhydroBody(i)->Pos[1]/PhydroBody(i)->Pos[0]);
                if(PhydroBody(i)->Pos[0] < 0) theta += M_PI;
                if(theta < 0) theta += 2*M_PI;

                PhydroBody(i)->Vel[0] += -v_tan*sin(theta);
                PhydroBody(i)->Vel[1] += +v_tan*cos(theta);
                if(GCHaloSpinVectorChange == 1){
                    double SpinTheta = GCHaloSpinTheta*M_PI/180.0;
                    double SpinPhi   = GCHaloSpinPhi*M_PI/180.0;
                    // for pos
                    // turn z-axis to x-axis
                    double Pos1[2] = {PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[2]};
                    PhydroBody(i)->Pos[0] = +Pos1[0]*cos(SpinTheta)+Pos1[1]*sin(SpinTheta);
                    PhydroBody(i)->Pos[2] = -Pos1[0]*sin(SpinTheta)+Pos1[1]*cos(SpinTheta);
                    // turn x-axis to y-axis
                    double Pos2[2] = {PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1]};
                    PhydroBody(i)->Pos[0] = +Pos2[0]*cos(SpinPhi)-Pos2[1]*sin(SpinPhi);
                    PhydroBody(i)->Pos[1] = +Pos2[0]*sin(SpinPhi)+Pos2[1]*cos(SpinPhi);
                    for(int k=0;k<3;k++)
                        PhydroBody(i)->PosP[k] = PhydroBody(i)->Pos[k];

                    // for vel
                    // turn z-axis to x-axis
                    double Vel1[2] = {PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[2]};
                    PhydroBody(i)->Vel[0] = +Vel1[0]*cos(SpinTheta)+Vel1[1]*sin(SpinTheta);
                    PhydroBody(i)->Vel[2] = -Vel1[0]*sin(SpinTheta)+Vel1[1]*cos(SpinTheta);
                    // turn x-axis to y-axis
                    double Vel2[2] = {PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1]};
                    PhydroBody(i)->Vel[0] = +Vel2[0]*cos(SpinPhi)-Vel2[1]*sin(SpinPhi);
                    PhydroBody(i)->Vel[1] = +Vel2[0]*sin(SpinPhi)+Vel2[1]*cos(SpinPhi);
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    for(int i=0;i<3;i++){
        GCInitPosForRemove[i] = InitPos[i];
    }
    //ParticleRemover(GCParticleRemoveCondition);
    ParticleRemoverArbitraryCriteria(GCParticleRemoveCondition);


    ActivateAllparticles();

    // UpdateTotalNumber();
    // ReConnectPointers();


    //write cloud.
    FILE *fp;
    char fcloud[MaxCharactersInLine];
    Snprintf(fcloud,"_cloud.%03d",MyID);
    FileOpen(fp,fcloud,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Tag == 0){
            fprintf(fp,"%d %g %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                    PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                    Phydro[i]->Tag);
        }
    }
    fclose(fp);
    fflush(NULL);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat _cloud.* > cloud.dat");
        fflush(NULL);
        system("rm -rf ./_cloud.*");
        fflush(NULL);
    }
    fflush(NULL);

    //write halo.
    Snprintf(fcloud,"_halo.%03d",MyID);
    FileOpen(fp,fcloud,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Tag == 1){
            fprintf(fp,"%d %g %g %g %g %g %g %d\n",PhydroBody(i)->GlobalID,
                    PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                    PhydroBody(i)->Vel[0],PhydroBody(i)->Vel[1],PhydroBody(i)->Vel[2],
                    Phydro[i]->Tag);
        }
    }
    fclose(fp);
    fflush(NULL);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat _halo.* > halo.dat");
        fflush(NULL);
        system("rm -rf ./_halo.*");
        fflush(NULL);
    }
    fflush(NULL);

    MPI_Barrier(MPI_COMM_WORLD);



    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"%s/Init.%02d.%02d",GCOutputDir,MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],Pbody[i]->Mass);
    }
    fclose(fp_init);
    {
        FILE *fp_runreport;
        char FileName[MaxCharactersInLine];
        Snprintf(FileName,"%s/InitReport.txt",GCOutputDir);
        FileOpen(fp_runreport,FileName,"w");
        fprintf(fp_runreport,"== Cloud Parameter ==\n");
        fprintf(fp_runreport," Particle number = %d\n",NumberCloud_t);
        fprintf(fp_runreport," Particle Mass  = %g [Earth Mass]\n",mgas*Pall.UnitMass/MEARTH_CGS);
        fprintf(fp_runreport," Softening Length  = %g [A.U.]\n",eps*Pall.UnitLength/AU_CGS);
        fprintf(fp_runreport,"\n\n");
        fprintf(fp_runreport,"== Halo Parameter ==\n");
        fprintf(fp_runreport," Particle number = %d\n",NumberHalo_t);
        fprintf(fp_runreport," Particle Mass  = %g [Earth Mass]\n",mhalo*Pall.UnitMass/MEARTH_CGS);
        fprintf(fp_runreport," Softening Length  = %g [A.U.]\n",epshalo*Pall.UnitLength/AU_CGS);
        fprintf(fp_runreport," Total Halo mass  = %g [Earth Mass]\n",GCReturnHaloMass(GCHaloTruncR));
        fprintf(fp_runreport,"\n\n");
        fprintf(fp_runreport,"== Sgr* Parameter ==\n");
        fprintf(fp_runreport," Mass  = %g [Msun]\n",GCMBH*Pall.UnitMass/MSUN_CGS);
        fprintf(fp_runreport," Softening Length  = %g [A.U.]\n",GCEPSBH*Pall.UnitLength/AU_CGS);
        fprintf(fp_runreport," Eat halo gas = %d\n",GCEATHALOGAS);

        fclose(fp_runreport);
    }
    fflush(NULL);

    GCWriteHaloProfile(Number,NumberCloud,NumberHalo);

    // ReConnectPointers();
    // UpdateTotalNumber();
    // UpdateTotalActiveNumber();

    dprintlmpi(NumberCloud_t);
    dprintlmpi(NumberHalo_t);

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    //Pall.Ns = 32;
    //Pall.Npm = 2;
    Pall.Ns = 64;
    Pall.Npm = 4;
    // Pall.Ns = 128;
    // Pall.Npm = 8;

    Pall.TEnd = (GCDay_end-GCDay_start)*YEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.001*Pall.TEnd;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
                Pall.TEnd,Pall.TEnd*Pall.UnitTime);
        fprintf(stderr,"The out put interval is %g in simulation unit, %g [s]\n",
                Pall.OutPutInterval,Pall.OutPutInterval*Pall.UnitTime);
    }


    if(MPIGetMyID() == MPI_ROOT_RANK){
        if(!CheckDir("./data"))
            MakeDir("./data");
    }

    strcpy(Pall.ASCIIFileName,"./data/GC.ASCII");
    strcpy(Pall.BaseFileName,"./data/GC");
    strcpy(Pall.RestartFileName,"./data/GC.dump");

    return;
}
#endif //TASK_GALACTIC_CENTER


void OutPutNFW(void){

    for(int k=0;k<MPIGetNumProcs();k++){
        if(k == MPIGetMyID()){
            double Temperature;
            FILE *fp_hydro,*fp_star;
            char fname_hydro[MaxCharactersInLine],fname_star[MaxCharactersInLine];
            sprintf(fname_hydro,"%s.Hydro.%02d.%02d",Pall.ASCIIFileName,MPIGetMyID(),MPIGetNumProcs());
            sprintf(fname_star,"%s.Star.%02d.%02d",Pall.ASCIIFileName,MPIGetMyID(),MPIGetNumProcs());
            FileOpen(fp_hydro,fname_hydro,"w");
            FileOpen(fp_star,fname_star,"w");

            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Type == TypeHydro){
                    fprintf(fp_hydro,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass,Pall.ConvertNumberDensityToCGS*PbodyHydro(i)->Rho,
                            PbodyHydroU(i)*Pall.ConvertUtoT,PbodyHydroKernel(i),
                            PbodyHydroZ(i),PbodyHydroZII(i),PbodyHydroZIa(i));

                } else if(Pbody[i]->Type == TypeStar){
                    fprintf(fp_star,"%ld %g %g %g %g %g %g %g %g %g %g %d %d\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass,
                            PbodyStarZ(i),PbodyStarZII(i),PbodyStarZIa(i),
                            PbodyStarTypeII(i),PbodyStarTypeIa(i));
                }
            }

            fclose(fp_hydro);
            fclose(fp_star);

            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

	return;
}

void OutPutAllParticlesInASCIIFormat(void){

    static double LastASCIIDataDump = 0;
    double CurrentTime = GetElapsedTime();
    // if(CurrentTime-LastASCIIDataDump < ASCIIDATA_DUMP_INTERVAL*60.e0 )
        // return;

    int write_flag = 0;
    if(CurrentTime-LastASCIIDataDump > ASCIIDATA_DUMP_INTERVAL*60.e0)
        write_flag ++;
    MPI_Allreduce(MPI_IN_PLACE,&write_flag,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(2*write_flag < MPIGetNumProcs()){
        return;
    }


    for(int k=0;k<MPIGetNumProcs();k++){
        if(k == MPIGetMyID()){
            double Temperature;
            FILE *fp_hydro,*fp_star,*fp_dm,*fp_sink;
            char fname_hydro[MaxCharactersInLine],fname_star[MaxCharactersInLine],
                 fname_dm[MaxCharactersInLine],fname_sink[MaxCharactersInLine];

            sprintf(fname_hydro,"%s.Hydro.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());
            sprintf(fname_star,"%s.Star.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());
            sprintf(fname_dm,"%s.DM.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());
            sprintf(fname_sink,"%s.Sink.%02d.%02d",Pall.ASCIIFileName,MPIGetNumProcs(),MPIGetMyID());

            if(Pall.Nhydro_t>0){
                FileOpen(fp_hydro,fname_hydro,"w");
            }
            if(Pall.Nstars_t>0){
                FileOpen(fp_star,fname_star,"w");
            }
            if(Pall.NDM_t>0){
                FileOpen(fp_dm,fname_dm,"w");
            }
            if(Pall.Nsink_t>0){
                FileOpen(fp_sink,fname_sink,"w");
            }

            for(int i=0;i<Pall.Ntotal;i++){
                if(Pbody[i]->Type == TypeHydro){
                    if(Pall.Nhydro_t>0)
                    fprintf(fp_hydro,"%ld %g %g %g %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass,Pall.ConvertNumberDensityToCGS*PbodyHydro(i)->Rho,
                            PbodyHydroU(i)*Pall.ConvertUtoT,PbodyHydroKernel(i),
                            PbodyHydroZ(i),PbodyHydroZII(i),PbodyHydroZIa(i));
                } else if(Pbody[i]->Type == TypeStar){
                    if(Pall.Nstars_t>0)
                    fprintf(fp_star,"%ld %g %g %g %g %g %g %g %g %g %g %d %d\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass,
                            PbodyStarZ(i),PbodyStarZII(i),PbodyStarZIa(i),
                            PbodyStarTypeII(i),PbodyStarTypeIa(i));
                } else if(Pbody[i]->Type == TypeSink){
                    if(Pall.Nsink_t>0)
                    fprintf(fp_sink,"%ld %g %g %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass,
                            PbodySink(i)->FormationTime,
                            PbodySink(i)->AM[0],PbodySink(i)->AM[1],PbodySink(i)->AM[2]);
                } else if(Pbody[i]->Type == TypeDM){
                    if(Pall.NDM_t>0)
                    fprintf(fp_dm,"%ld %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                            Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                            Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                            Pbody[i]->Mass);
                }
            }

            if(Pall.Nhydro_t>0){
                fclose(fp_hydro);
            }
            if(Pall.Nstars_t>0){
                fclose(fp_star);
            }
            if(Pall.NDM_t>0){
                fclose(fp_dm);
            }
            if(Pall.Nsink_t>0){
                fclose(fp_sink);
            }
            fflush(NULL);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    LastASCIIDataDump = CurrentTime;

	return;
}

void ParticleDistributionGenetator(int number, unsigned long int seed){

    int MyID = MPIGetMyID();
    int AllocationSize = number;

    memset(&Pall,0,sizeof(struct StructPall));
    InitializeRandomGenerator(1977+MyID);
    //Pall.Lbox = 1.e0;
    //Pall.Lboxh = Pall.Lbox/2.e0;

    Pall.Ntotal = number;
    PbodySize = AllocationSize;

    PbodyElements = malloc(AllocationSize*sizeof(StructPbody));
    Pbody = malloc(PbodySize*sizeof(StructPbodyptr));
    memset(PbodyElements,0,AllocationSize*sizeof(StructPbody));

    for(int i=0;i<AllocationSize-1;i++)
        PbodyElements[i].Next = &(PbodyElements[i+1]);

    for(int i=0;i<number;i++){
        Pbody[i] = PbodyElements+i;
        PbodyElements[i].Active = ON;
        PbodyElements[i].Use = ON;
        PbodyElements[i].Type = TypeDM;
        PbodyElements[i].GlobalID = i+number*MyID;
        PbodyElements[i].Pos[0] = 0.5*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        PbodyElements[i].Pos[1] = 0.5*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        PbodyElements[i].Pos[2] = 0.5*(2.0*gsl_rng_uniform(RandomGenerator)-1.0);
        PbodyElements[i].Mass = 1.e0;
        PbodyElements[i].Eps = 1.e-5;
        PbodyElements[i].Baryon = NULL;
        Pall.NDM ++;
    }

    if(MPIGetMyID() == 0){
        for(int i=0;i<number;i+=500){
            fprintf(stderr,"%d %ld %g %d\n",
                    i,Pbody[i]->GlobalID,Pbody[i]->Pos[0],Pbody[i]->Type);
        }
    }

    /** Initializetions **/
    Pall.AdaptiveSofteningFactor = 1.0;
    /** Initializetions **/
    
    return;
}

#ifdef TASK_TEST_STROMGRENSPHERE
static double LyZSolar = 0.243E+47; // Photon number per solar mass of an SSP particle.
void InitTestStromgrenSphere(const double nH, const double Radius, const int Nhydro, const int Nstars, double Pos[restrict][3], double LyAlpha[restrict]){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitLength = PC_CGS;
    Pall.UnitTime = MEGAYEAR_CGS;
    Pall.UnitMass = MSUN_CGS;

    Pall.GravConst = GetUnitGravitationalConstant();

    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();

    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    // 
    int m = (int)log2(Nhydro);
    int Nhydro_p2 = 1<<m;
    dprintlmpi(Nhydro_p2);

    double ParticleMass = (4*M_PI/3.0)*(nH*PROTON_MASS_CGS/0.76)*CUBE(Radius*PC_CGS)
                            /MSUN_CGS/(double)Nhydro_p2;
    fprintf(stderr,"ParticleMass = %g, Total Mass = %g\n",ParticleMass,ParticleMass*Nhydro_p2);


    if(MPIGetMyID() == MPI_ROOT_RANK){
        Pall.Nhydro = Nhydro_p2/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nstars = Pall.Nstars_t = Nstars;
        Pall.Nsink = Pall.Nsink_t = 0;
        Pall.Ntotal = Pall.Nhydro + Pall.Nstars;
        Pall.Ntotal_t = Nhydro_p2 + Nstars;
    } else {
        Pall.Nhydro = Nhydro_p2/NProcs;
        Pall.Nhydro_t = Pall.Nhydro*NProcs;
        Pall.Nstars = 0;
        Pall.Nstars_t = Nstars;
        Pall.Nsink = Pall.Nsink_t = 0;
        Pall.Ntotal = Pall.Nhydro + Pall.Nstars;
        Pall.Ntotal_t = Nhydro_p2 + Nstars;
    }
    fprintf(stderr,"[%02d] %ld %ld %ld %ld %ld | %ld %ld %ld %ld %ld\n",MPIGetMyID(),
            Pall.Nhydro,Pall.Nstars,Pall.Nsink,Pall.NDM,Pall.Ntotal,
            Pall.Nhydro_t,Pall.Nstars_t,Pall.Nsink_t,Pall.NDM_t,Pall.Ntotal_t);


    int AllocationSize = Pall.Ntotal; 
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(Pall.Nhydro);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        GenerateStructPstar(Nstars);
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }
    if(MPIGetMyID() == MPI_ROOT_RANK){
        for(int i=0;i<Nstars;i++){
            Pbody[Pall.Nhydro+i]->Baryon = (void *)(Pstar[i]);
            Pstar[i]->Body = Pbody[Pall.Nhydro+i];
        }
    }

    double eps = 0.1*PC_CGS/Pall.UnitLength;
    double Uinit = 100*Pall.ConvertTtoU;
    double cs = sqrt(Pall.GGm1*Uinit);

#if 1
    int count = 0;
    int count_passall = 0;
    while(count_passall < Pall.Nhydro_t){
        double Pos[] ={Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0),
                       Radius*(2.0*gsl_rng_uniform(RandomGenerator)-1.0)};
        // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
        double r = NORM(Pos);
        if(r<Radius){
            // fprintf(stderr,"%g %g %g\n",Pos[0],Pos[1],Pos[2]);
            if(count_passall%NProcs == MyID){
                Pbody[count]->Active = ON;
                Pbody[count]->Use = ON;
                Pbody[count]->Type = TypeHydro;
                Pbody[count]->GlobalID = count_passall;

                Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
                Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
                Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Pos[2];

                Pbody[count]->Vel[0] = Pbody[count]->Vel[1] = Pbody[count]->Vel[2] = 1.e-30;
                Pbody[count]->Acc[0] = Pbody[count]->Acc[1] = Pbody[count]->Acc[2] = 0.e0;

                Pbody[count]->Eps = eps;
                Pbody[count]->Mass = ParticleMass;

                PbodyHydro(count)->Mass = Pbody[count]->Mass = ParticleMass;
                PbodyHydro(count)->Z = 0.02;
                PbodyHydro(count)->ZII = 0.02;
                PbodyHydro(count)->ZIa = 0.02;
#if (UseSFModelSpawn) 
                PbodyHydro(count)->SpawnMass = ParticleMass/((double)MaxSpawnTimes);
#endif

#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
                PbodyHydro(count)->G0 = 10;
                PbodyHydro(count)->fH2 = 0.4;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

                PbodyHydro(count)->DQheat = 0.e0;
                PbodyHydro(count)->Use = ON;
                PbodyHydro(count)->Kernel = eps;
                PbodyHydro(count)->U = Uinit;

                count ++;
            }
            count_passall ++;
        }
    }


#else
    int count = 0;
    for(int i=0;i<Pall.Nhydro;i++){
        double R = Radius*gsl_rng_uniform(RandomGenerator);
        double costh = 2.0*gsl_rng_uniform(RandomGenerator)-1;
        double phi = 360*gsl_rng_uniform(RandomGenerator);
        double phi_rad = (M_PI/180.0)*phi;
        if(i%NProcs == MyID){
            // fprintf(stderr,"%d\n",i);
            double Pos[3] = {R*sqrt(1-SQ(costh))*cos(phi_rad),
                            R*sqrt(1-SQ(costh))*sin(phi_rad),
                            R*costh};
            
            Pbody[count]->Active = ON;
            Pbody[count]->Use = ON;
            Pbody[count]->Type = TypeHydro;
            Pbody[count]->GlobalID = i;

            Pbody[count]->PosP[0] = Pbody[count]->Pos[0] = Pos[0];
            Pbody[count]->PosP[1] = Pbody[count]->Pos[1] = Pos[1];
            Pbody[count]->PosP[2] = Pbody[count]->Pos[2] = Pos[2];

            Pbody[count]->Vel[0] = Pbody[count]->Vel[1] = Pbody[count]->Vel[2] = 1.e-30;
            Pbody[count]->Acc[0] = Pbody[count]->Acc[1] = Pbody[count]->Acc[2] = 0.e0;

            Pbody[count]->Eps = eps;
            Pbody[count]->Mass = ParticleMass;

            PbodyHydro(count)->Mass = Pbody[count]->Mass = ParticleMass;
            PbodyHydro(count)->Z = 0.02;
            PbodyHydro(count)->ZII = 0.02;
            PbodyHydro(count)->ZIa = 0.02;
#if (UseSFModelSpawn) 
            PbodyHydro(count)->SpawnMass = ParticleMass/((double)MaxSpawnTimes);
#endif

#ifdef USE_SPAANS2008_COOLING_FUNCTIONS
            PbodyHydro(count)->G0 = 10;
            PbodyHydro(count)->fH2 = 0.4;
#endif //USE_SPAANS2008_COOLING_FUNCTIONS

            PbodyHydro(count)->DQheat = 0.e0;
            PbodyHydro(count)->Use = ON;
            PbodyHydro(count)->Kernel = eps;
            PbodyHydro(count)->U = Uinit;

            count ++;
        }
    }
#endif
    // dprintlmpi(count);
    // dlprintlmpi(Pall.Nhydro);

#if 0
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }
#endif

// star particles.
    if(MPIGetMyID() == MPI_ROOT_RANK){
        for(int i=0;i<Pall.Nstars;i++){
            int index = Pall.Nhydro+i;
            Pbody[index]->Active = ON;
            Pbody[index]->Use = ON;
            Pbody[index]->Type = TypeStar;
            Pbody[index]->GlobalID = Pall.Nhydro_t+i;

            Pbody[index]->Pos[0] = Pos[i][0];
            Pbody[index]->Pos[1] = Pos[i][1];
            Pbody[index]->Pos[2] = Pos[i][2];
            Pbody[index]->Vel[0] = Pbody[index]->Vel[1] = Pbody[index]->Vel[2] = 0.e0;
            Pbody[index]->Acc[0] = Pbody[index]->Acc[1] = Pbody[index]->Acc[2] = 0.e0;
            Pbody[index]->Eps = eps;
            Pbody[index]->Mass = LyAlpha[i]/LyZSolar;

            PbodyStar(index)->Use = ON;
            PbodyStar(index)->ParentGlobalID = Pbody[index]->GlobalID;
            PbodyStar(index)->FormationTime = 0.e0;
#ifdef PRESERVE_SNII_EVENTRATE //{
            PbodyStar(index)->TypeIIProb = true; 
#endif // PRESERVE_SNII_EVENTRATE //}
            PbodyStar(index)->InitialMass = Pbody[index]->Mass; 
            fprintf(stderr,"Stellar Mass is %g [Msun]\n",PbodyStar(index)->InitialMass*Pall.UnitMass/MSUN_CGS);
        }
    }

    ActivateAllparticles();


    FILE *fp_init;
    char Fname[MaxCharactersInLine];
    Snprintf(Fname,"./init.%02d.%02d",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp_init,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp_init,"%ld %g %g %g %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2],
                Pbody[i]->Acc[0],Pbody[i]->Acc[1],Pbody[i]->Acc[2]);
    }
    fclose(fp_init);

#if 1
    int Nbin = 30;
    double Rad[Nbin];
    double Mass[Nbin];
    int Pnum[Nbin];
    double dr = Radius/(double)Nbin;
    gprintlmpi(dr);
    for(int i=0;i<Nbin;i++){
        Rad[i] = (i+0.5)*dr;
        Mass[i] = 0.e0;
        Pnum[i] = 0;
    }

    for(int i=0;i<Pall.Nhydro;i++){
        double r = NORM(PhydroBody(i)->Pos);
        int int_r = r/dr;
        if(int_r < 0) continue;
        if(int_r >= Nbin) continue;

        Mass[int_r] += PhydroBody(i)->Mass;
        Pnum[int_r] ++;
    }
    double GlobalMass[Nbin];
    MPI_Allreduce(Mass,GlobalMass,Nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    int GlobalPnum[Nbin];
    MPI_Allreduce(Pnum,GlobalPnum,Nbin,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(MyID == MPI_ROOT_RANK){
        FILE *fp;
        FileOpen(fp,"./rad.dist","w");

        double m = 0.e0; 
        // output mean (number) density.
        for(int i=0;i<Nbin;i++){
            m += Mass[i];
            double density = m*Pall.UnitMass/((4*M_PI/3.0)*CUBE(dr*(i+1)*PC_CGS));
            double Ndensity = 0.76*density/PROTON_MASS_CGS;
            // fprintf(stderr,"Mass %g, R %g, density %g, nH %g\n",
                    // m,dr*(i+1),density,Ndensity);
            fprintf(fp,"%g %g %d %g %g\n",Rad[i],GlobalMass[i],GlobalPnum[i],density,Ndensity);
        }
        fclose(fp);
    }
#endif

    Pall.NActivesHydro = Pall.Nhydro;
    Pall.NActivesStars = Pall.Nstars;
    Pall.NActivesSink = Pall.Nsink;
    Pall.NActivesDM = Pall.NDM;
    Pall.NActives = Pall.NActivesHydro+Pall.NActivesStars+Pall.NActivesSink+Pall.NActivesDM;

    ReConnectPointers();
    UpdateTotalNumber();
    UpdateTotalActiveNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;

#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.TEnd = 100*MEGAYEAR_CGS/Pall.UnitTime;
    Pall.TCurrent = 0.e0;
    Pall.Redshift = Pall.InitialRedshift = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    if(MPIGetMyID() == MPI_ROOT_RANK){
        fprintf(stderr,"The end time of this simulation is %g in simulation unit, %g [s]\n",
                Pall.TEnd,Pall.TEnd*Pall.UnitTime);
    }


    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.0001*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/HII.ASCII");
    strcpy(Pall.BaseFileName,"./data/HII");
    strcpy(Pall.RestartFileName,"./data/HII.dump");

    return;
}

#endif //TASK_TEST_STROMGRENSPHERE

#ifdef TASK_HYDROSTATIC
static void WriteRunLog(const int FlagEqualMass, const double DensityRatio){

    FILE *fp;
    char fname[MaxCharactersInLine];

    MakeDir("./data");

    sprintf(fname,"./data/hs.log");
    FileOpen(fp,fname,"w");

    fprintf(fp,"%d\n",FlagEqualMass);
    fprintf(fp,"%g\n",DensityRatio);

    fclose(fp);

    return ;
}

/* 
 * This is the initial particle distribution generator for the hydrostatic test.
 * The reference is Saitoh & Makino (2013)
 */
void InitHydroStatic(const int NGrid, const int FlagEqualMass, const double DensityRatio){

    WriteRunLog(FlagEqualMass,DensityRatio);

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);


    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.GravConst = 
    Pall.ConvertTtoU = Pall.ConvertUtoT = 
    Pall.ConvertDensityToCGS = Pall.ConvertNumberDensityToCGS = 1.0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else // USE_VARIABLE_ALPHA //}//{
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    int Base = (int)(sqrt(DensityRatio)+0.5);
    dprintl(Base);

    double rho_h = SQ(Base);
    double rho_l = 1.0;
    double p_h = 2.5;
    double p_l = 2.5;
    double u_h = p_h/(Pall.Gm1*rho_h);
    double u_l = p_l/(Pall.Gm1*rho_l);

    double dx = 1.0/(double)NGrid;
    int AllocationSize = 0; int _counter = 0;
    for(int i=0;i<NGrid;i++){
        for(int k=0;k<NGrid;k++){
            if(_counter%NProcs == MyID){
                double x = dx*(i+0.5);
                double y = dx*(k+0.5);
                if((x>0.25)&&(x<0.75)&&(y>0.25)&&(y<0.75)){
                    if(FlagEqualMass == 1){
                        AllocationSize += SQ(Base);
                    } else {
                        AllocationSize ++;
                    }
                }  else {
                    AllocationSize ++;
                }
            }
            _counter ++;
        }
    }
    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double mgas = 1.0/(SQ(NGrid));

    int mycount = 0; int counter = 0;
    for(int i=0;i<NGrid;i++){
        for(int k=0;k<NGrid;k++){
            double x = dx*(i+0.5);
            double y = dx*(k+0.5);
            if((x>0.25)&&(x<0.75)&&(y>0.25)&&(y<0.75)){
                if(FlagEqualMass){
                    double dxx = dx/Base;
                    for(int ii=0;ii<Base;ii++){
                        for(int kk=0;kk<Base;kk++){
                            double xx = x + dxx*(ii-0.5);
                            double yy = y + dxx*(kk-0.5);

                            if(counter%NProcs == MyID){
                                PhydroBody(mycount)->Active = ON;
                                PhydroBody(mycount)->Use = ON;
                                PhydroBody(mycount)->Type = TypeHydro;
                                PhydroBody(mycount)->GlobalID = counter;


                                PhydroBody(mycount)->Pos[0] = xx;
                                PhydroBody(mycount)->Pos[1] = yy;
                                PhydroBody(mycount)->Pos[2] = 0.5;
                                PhydroBody(mycount)->Vel[0] =
                                PhydroBody(mycount)->Vel[1] =
                                PhydroBody(mycount)->Vel[2] = 0.e0;
                                PhydroBody(mycount)->Mass = 
                                    Phydro[mycount]->Mass = mgas;

                                Phydro[mycount]->Rho = rho_h;
                                Phydro[mycount]->Kernel = 2*dxx;
                                Phydro[mycount]->U = u_h;

                                Phydro[mycount]->Du =
                                Phydro[mycount]->HydroAcc[0] =
                                Phydro[mycount]->HydroAcc[1] =
                                Phydro[mycount]->HydroAcc[2] = 0.e0;

                                Phydro[mycount]->Use = ON; 
                                PhydroBody(mycount)->Eps = 1.0;
#ifdef USE_PARTICLE_TAG
                                Phydro[mycount]->Tag = 1; 
#endif // USE_PARTICLE_TAG
                                mycount ++;
                            }
                            counter ++;
                        }
                    }
                } else {
                    if(counter%NProcs == MyID){
                        PhydroBody(mycount)->Active = ON;
                        PhydroBody(mycount)->Use = ON;
                        PhydroBody(mycount)->Type = TypeHydro;
                        PhydroBody(mycount)->GlobalID = counter;

                        PhydroBody(mycount)->Pos[0] = x;
                        PhydroBody(mycount)->Pos[1] = y;
                        PhydroBody(mycount)->Pos[2] = 0.5;
                        PhydroBody(mycount)->Vel[0] =
                        PhydroBody(mycount)->Vel[1] =
                        PhydroBody(mycount)->Vel[2] = 0.e0;
                        PhydroBody(mycount)->Mass = 
                            Phydro[mycount]->Mass = mgas*SQ(Base);

                        Phydro[mycount]->Rho = rho_h;
                        Phydro[mycount]->Kernel = 2*dx;
                        Phydro[mycount]->U = u_h;

                        Phydro[mycount]->Du =
                        Phydro[mycount]->HydroAcc[0] =
                        Phydro[mycount]->HydroAcc[1] =
                        Phydro[mycount]->HydroAcc[2] = 0.e0;

                        Phydro[mycount]->Use = ON; 
                        PhydroBody(mycount)->Eps = 1.0;
#ifdef USE_PARTICLE_TAG
                        Phydro[mycount]->Tag = 1; 
#endif // USE_PARTICLE_TAG
                        mycount ++;
                    }
                    counter ++;
                }
            }  else {
                if(counter%NProcs == MyID){
                    PhydroBody(mycount)->Active = ON;
                    PhydroBody(mycount)->Use = ON;
                    PhydroBody(mycount)->Type = TypeHydro;
                    PhydroBody(mycount)->GlobalID = counter;

                    PhydroBody(mycount)->Pos[0] = x;
                    PhydroBody(mycount)->Pos[1] = y;
                    PhydroBody(mycount)->Pos[2] = 0.5;
                    PhydroBody(mycount)->Vel[0] =
                    PhydroBody(mycount)->Vel[1] =
                    PhydroBody(mycount)->Vel[2] = 0.e0;
                    PhydroBody(mycount)->Mass = 
                        Phydro[mycount]->Mass = mgas;

                    Phydro[mycount]->Rho = rho_l;
                    Phydro[mycount]->Kernel = 2*dx;
                    Phydro[mycount]->U = u_l;

                    Phydro[mycount]->Du =
                    Phydro[mycount]->HydroAcc[0] =
                    Phydro[mycount]->HydroAcc[1] =
                    Phydro[mycount]->HydroAcc[2] = 0.e0;

                    Phydro[mycount]->Use = ON; 
                    PhydroBody(mycount)->Eps = 1.0;
#ifdef USE_PARTICLE_TAG
                    Phydro[mycount]->Tag = 0; 
#endif // USE_PARTICLE_TAG
                    mycount ++;
                }
                counter ++;
            }
        }
    }

    dprintl(mycount);
    Pall.Ntotal = Pall.Nhydro = mycount;
    MPI_Allreduce(&Pall.Ntotal,&Pall.Ntotal_t,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&Pall.Nhydro,&Pall.Nhydro_t,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);

    ActivateAllparticles();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    //Pall.Ns = 64;
    //Pall.Ns = 128;
    Pall.Npm = 2;

    Pall.TEnd = 8.0;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/80;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/hs.ASCII");
    strcpy(Pall.BaseFileName,"./data/hs");
    strcpy(Pall.RestartFileName,"./data/hs.dump");

    return ;
}

#endif //TASK_HYDROSTATIC

#ifdef TASK_KELVINHELMHOLTZ_INSTABILITY //{
static void WriteRunLog(const int NGrid, const int mode, const int seed,
        const double density_contrast, const double Tkh, const double Tend){

    FILE *fp;
    char fname[MaxCharactersInLine];
    sprintf(fname,"./data/kh.log");
    FileOpen(fp,fname,"w");

    fprintf(fp,"%d\n",NGrid);
    fprintf(fp,"%d\n",mode);
    fprintf(fp,"%d\n",seed);
    fprintf(fp,"%g\n",density_contrast);
    fprintf(fp,"%g\n",Tkh);
    fprintf(fp,"%g\n",Tend);

    fclose(fp);

    return ;
}

/* 
 * The kelvin-Helmoltz instability tests.
 * References are Price, D. (2008), Read, J. et al. (2010), and Saitoh & Makino
 * (2013).
*/
void InitKH(const int NGrid, const int mode, const int seed, 
        const double density_contrast, const double TEnd){
    // mode: 0->Price type, 1->Read type initial condition.
    // seed: 0->Price tyep, 1->Read type initial velocity perturbation.
    // density_constrast is the density contrast between two phase.
    // TEnd is the time to finish this run.

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977);

    Pall.UnitTime = 1;
    Pall.UnitLength = 1;
    Pall.UnitMass = 1;
    Pall.Lbox[0] = Pall.Lbox[1] = Pall.Lbox[2] = 1.e0;
    Pall.Lboxh[0] = Pall.Lboxh[1] = Pall.Lboxh[2] = 0.5;

    Pall.GravConst = 
    Pall.ConvertTtoU = Pall.ConvertUtoT = 
    Pall.ConvertDensityToCGS = Pall.ConvertNumberDensityToCGS = 1.0;

    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    double rho_h = density_contrast;
    gprintlmpi(density_contrast);
    gprintlmpi(TEnd);
    double rho_l = 1.0;
    double p_h = 2.5;
    double p_l = 2.5;
    double v_h = 0.5;
    double v_l = -0.5;
    double u_h = p_h/(Pall.Gm1*rho_h);
    double u_l = p_l/(Pall.Gm1*rho_l);

    double cs_h = sqrt(Pall.GGm1*u_h);
    double cs_l = sqrt(Pall.GGm1*u_l);

    fprintf(stderr,"cs_h = %g, cs_l = %g\n",cs_h,cs_l);

    int AllocationSize = 0;
    int ngrid_h,ngrid_l,n_h,n_l;
    double mgas_h,mgas_l;
    if(mode == 0){ // equal mass.
        ngrid_h = NGrid;
        ngrid_l = sqrt(SQ(ngrid_h) *rho_l/rho_h);
        n_h = 0.5 * SQ(ngrid_h);
        n_l = 0.5 * SQ(ngrid_l);
        for(int i=0;i<n_h;i++){
            if(i%NProcs == MyID){
                AllocationSize ++;
            }
        }
        for(int i=0;i<n_l;i++){
            if(i%NProcs == MyID){
                AllocationSize ++;
            }
        }
        mgas_h = mgas_l = 0.5 * (rho_h + rho_l) / ((double)(n_h + n_l));
    } else if(mode == 1){ // equal sep.
        ngrid_h = NGrid;
        ngrid_l = NGrid;
        n_h = 0.5 * SQ(ngrid_h);
        n_l = 0.5 * SQ(ngrid_l);
        for(int i=0;i<n_h;i++){
            if(i%NProcs == MyID){
                AllocationSize ++;
            }
        }
        for(int i=0;i<n_l;i++){
            if(i%NProcs == MyID){
                AllocationSize ++;
            }
        }
        mgas_h = 0.5 * (rho_h) / ((double)n_h);
        mgas_l = 0.5 * (rho_l) / ((double)n_l);
    } else {
        fprintf(stderr,"Mode setting error!\n");
        fflush(NULL);
        exit(EXIT_SUCCESS);
    }

    double KH_Width,KH_Lambda,KH_Amp;
    if(seed == 0){
        KH_Width  = 0.025;
        KH_Lambda = 1.0/6.0;
        KH_Amp    = 0.025;
    } else if(seed == 1){
        KH_Width  =  0.025;
        KH_Lambda =  1.0/2.0;
        KH_Amp    =  1.0/8.0;
    }

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    int mycount = 0; int counter = 0;
    int half_grid = 0.5 * ngrid_h;
    double dx_h = sqrt(0.5/((double)n_h));
    double offset_h = 0.5*(0.5-dx_h*half_grid);
    for (int i = 0; i != ngrid_h; ++i) {
        for (int j = 0; j != half_grid; ++j) {
            if(counter%NProcs == MyID){
                PhydroBody(mycount)->Active = ON;
                PhydroBody(mycount)->Use = ON;
                PhydroBody(mycount)->Type = TypeHydro;
                PhydroBody(mycount)->GlobalID = counter;

                double x = (0.5 + i) * dx_h;
                //double y = 0.25 + (0.5 + j) * dx_h;
                double y = 0.25+offset_h + (0.5 + j) * dx_h;
                Phydro[mycount]->Kernel = 2*dx_h;
                PhydroBody(mycount)->Pos[0] = x;
                PhydroBody(mycount)->Pos[1] = y;
                PhydroBody(mycount)->Pos[2] = 0.5;

                if(seed == 0){
                    PhydroBody(mycount)->Vel[0] = v_h;
                    double vy = 0;
                    if (fabs(y - 0.5 - 0.25) < KH_Width)
                        vy = KH_Amp * sin(-2 * M_PI * (x+0.5)/KH_Lambda);
                    else if (fabs(y - 0.5 + 0.25) < KH_Width)
                        vy = KH_Amp * sin(2 * M_PI * (x+0.5)/KH_Lambda);
                    PhydroBody(mycount)->Vel[1] = vy;
                } else if(seed == 1){
                    PhydroBody(mycount)->Vel[0] = v_h;
                    PhydroBody(mycount)->Vel[1] = KH_Amp*
                        (sin(2.0*M_PI*(x+KH_Lambda/2.0)/KH_Lambda)*exp(-SQ(10*(y-0.75)))
                        -sin(2.0*M_PI*x/KH_Lambda)*exp(-SQ(10*(y-0.25))));
                }
                PhydroBody(mycount)->Vel[2] = 0.e0;

                PhydroBody(mycount)->Mass = Phydro[mycount]->Mass = mgas_h;
                Phydro[mycount]->U = u_h;
                Phydro[mycount]->Rho = rho_h;

                Phydro[mycount]->Du =
                Phydro[mycount]->HydroAcc[0] =
                Phydro[mycount]->HydroAcc[1] =
                Phydro[mycount]->HydroAcc[2] = 0.e0;

                Phydro[mycount]->Use = ON; 
                PhydroBody(mycount)->Eps = 1.0;
#ifdef USE_PARTICLE_TAG
                Phydro[mycount]->Tag = 1; 
#endif // USE_PARTICLE_TAG
                mycount ++;
            }
            counter ++;
        }
    }


    half_grid = 0.5 * ngrid_l;
    double dx_l = sqrt(0.5/((double)n_l));
    double offset_l = 0.5*(0.5-dx_l*half_grid);
    for (int i = 0; i != ngrid_l; ++i) {
        for (int j = 0; j != half_grid; ++j){
            if(counter%NProcs == MyID){
                PhydroBody(mycount)->Active = ON;
                PhydroBody(mycount)->Use = ON;
                PhydroBody(mycount)->Type = TypeHydro;
                PhydroBody(mycount)->GlobalID = counter;

                double x = (0.5 + i) * dx_l;
                //double y = (0.5 + j) * dx_l;
                //double y = (1.0 + j) * dx_l;
                // double y = -0.25 + (0.5 + j) * dx_l;
                double y = -0.25+offset_l + (0.5 + j) * dx_l;
                if (y > 0.25)
                    y += 0.5;
                if (y < 0.0)
                    y+= 1.0;
                Phydro[mycount]->Kernel = 2*dx_l;
                PhydroBody(mycount)->Pos[0] = x;
                PhydroBody(mycount)->Pos[1] = y;
                PhydroBody(mycount)->Pos[2] = 0.5;

                if(seed == 0){
                    PhydroVel(mycount)[0] = v_l;
                    double vy = 0;
                    if (fabs(y - 0.5 - 0.25) < KH_Width)
                        vy = KH_Amp * sin(-2 * M_PI * (x+0.5)/KH_Lambda);
                    else if (fabs(y - 0.5 + 0.25) < KH_Width)
                        vy = KH_Amp * sin(2 * M_PI * (x+0.5)/KH_Lambda);
                    PhydroVel(mycount)[1] = vy;
                } else if(seed ==1){
                    PhydroVel(mycount)[0] = v_l;
                    PhydroVel(mycount)[1] = KH_Amp*
                        (sin(2.0*M_PI*(x+KH_Lambda/2.0)/KH_Lambda)*exp(-SQ(10*(y-0.75)))
                        -sin(2.0*M_PI*x/KH_Lambda)*exp(-SQ(10*(y-0.25))));
                }
                PhydroVel(mycount)[2] = 0.e0;

                PhydroBody(mycount)->Mass = PhydroBody(mycount)->Mass = mgas_l;
                Phydro[mycount]->U = u_l;
                Phydro[mycount]->Rho = rho_l;

                Phydro[mycount]->Du =
                Phydro[mycount]->HydroAcc[0] =
                Phydro[mycount]->HydroAcc[1] =
                Phydro[counter]->HydroAcc[2] = 0.e0;

                Phydro[mycount]->Use = ON; 
                PhydroBody(mycount)->Eps = 1.0;
#ifdef USE_PARTICLE_TAG
                Phydro[mycount]->Tag = 0; 
#endif // USE_PARTICLE_TAG

                mycount ++;
            }
            counter ++;
        }
    }
    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = counter;


    ActivateAllparticles();

#if 0
    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
        Pbody[i]->Active = ON;
        Pbody[i]->Acc[0] = Pbody[i]->Acc[1] = Pbody[i]->Acc[2] = 0.e0;
    }
#endif

    double TNorm = ((rho_l+rho_h)*(KH_Lambda))/(sqrt(rho_l*rho_h)*(v_h-v_l));
    double _TEnd = TEnd*TNorm;
    if(MPIGetMyID() == MPI_ROOT_RANK)
        fprintf(stderr,"KH time scale is %g, Tend = %g\n",TNorm,TEnd);

    WriteRunLog(NGrid,mode,seed,density_contrast,TNorm,_TEnd);

    // ReConnectPointers();
    // UpdateTotalNumber();
    // UpdateTotalActiveNumber();

    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;


    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = _TEnd;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = Pall.TEnd/200;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/KH.ASCII");
    strcpy(Pall.BaseFileName,"./data/KH");
    strcpy(Pall.RestartFileName,"./data/KH.dump");


#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    return ;
}
#endif // TASK_KELVINHELMHOLTZ_INSTABILITY //}

#ifdef TASK_TEST_1D_THERMAL_CONDUCTIVITY //{
void Init1DThermalConductivity(const int Number, const double T_JUMP){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    Pall.UnitLength = 1;
    Pall.UnitTime = 1; // Velocity unit is set to 1km/s
    Pall.UnitMass = 1;
    Pall.TCMB = 0.0;

    Pall.Lbox[0] = 1.e0;
    Pall.Lbox[1] = 0.5;
    Pall.Lbox[2] = 0.5;
    Pall.Lboxh[0] = 0.5*Pall.Lbox[0];
    Pall.Lboxh[1] = 0.5*Pall.Lbox[1];
    Pall.Lboxh[2] = 0.5*Pall.Lbox[2];
    Pall.BoxCenter[0] = 0.0;
    Pall.BoxCenter[1] = 0.25;
    Pall.BoxCenter[2] = 0.25;

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else // USE_VARIABLE_ALPHA //}//{
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    // Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    // Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertTtoU = 1.e-3; 
    Pall.ConvertUtoT = 1.e+3; 
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    int NLeft = Number/2;
    int NRight = Number/2;
    int Ntotal = NLeft+NRight;

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Ntotal;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    dprintlmpi(AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }


    double dx = 0.5/((double)NLeft);
    double mass = 2.4e4/((double)Ntotal);

    dprintlmpi(Ntotal);
    int mycount = 0;
    for(int i=0;i<Ntotal;i++){
        if(i%NProcs == MyID){
            PhydroBody(mycount)->Active = ON;
            PhydroBody(mycount)->Use = ON;
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->GlobalID = i;
            if(i<NLeft){
                PhydroPos(mycount)[0] = -5.e-1+i*dx;
                Phydro[mycount]->Kernel = 2.0*dx;
                Phydro[mycount]->Rho = 2.4e+4;
                // Phydro[mycount]->U = 3.0e+15*50.e0/(Pall.Gm1*Phydro[mycount]->Rho);
                // Phydro[mycount]->U = 50.0;
            }else{
                PhydroPos(mycount)[0] = (i-NLeft)*dx;
                Phydro[mycount]->Kernel = 2.0*dx;
                Phydro[mycount]->Rho = 2.4e+4;
                //Phydro[mycount]->U = 3.0e+15*50.e0/(100)/(Pall.Gm1*Phydro[mycount]->Rho);
                // Phydro[mycount]->U = 50.0/(100);
            }

            double t = 0.000001;
            double Kappa = 1.0/(Phydro[mycount]->Rho*Pall.ConvertTtoU);
            double X = PhydroPos(mycount)[0]/2.0/sqrt(Kappa*t);
            double U = (50.0/2.0)*erfc(X)+(50.0/T_JUMP/2.0)*erfc(-X);
            Phydro[mycount]->U = U;
            PhydroPos(mycount)[0] += 0.5;


#ifdef EVALUATE_KERNEL_BY_ITERATION
            Phydro[mycount]->Kernel = 1.e0/((double)Number);
#endif
            PhydroBody(mycount)->Mass = PhydroMass(mycount) = mass;
            PhydroPos(mycount)[1] = PhydroPos(mycount)[2] = 0.e0;
            PhydroVel(mycount)[0] = PhydroVel(mycount)[1] = PhydroVel(mycount)[2] = 0.e0;
            Phydro[mycount]->HydroAcc[0] = Phydro[mycount]->HydroAcc[1] = Phydro[mycount]->HydroAcc[2] = 0.e0;
            Phydro[mycount]->Use = ON; 
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->Eps = 1.0;
            mycount ++;
        }
    }
    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = Ntotal;

    ActivateAllparticles();

    FILE *fp;
    char fcloud[MaxCharactersInLine];
    Snprintf(fcloud,"_hydro.%03d",MyID);
    FileOpen(fp,fcloud,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%d %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->U);
    }
    fclose(fp);
    fflush(NULL);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat _hydro.* > hydro.dat");
        fflush(NULL);
        system("rm -rf ./_hydro.*");
        fflush(NULL);
    }
    fflush(NULL);


    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 6;
    Pall.Npm = 2;

    Pall.TEnd = 0.024;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/TC1D.ASCII");
    strcpy(Pall.BaseFileName,"./data/TC1D");
    strcpy(Pall.RestartFileName,"./data/TC1D.dump");


    return;
}
#endif // TASK_TEST_1D_THERMAL_CONDUCTIVITY //}

#ifdef TASK_TEST_3D_THERMAL_CONDUCTIVITY //{
void Init3DThermalConductivity(const int Number, const int mode, const double T_JUMP){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));
    memset(&TimingResults,0,sizeof(struct StructTimingResults));
    InitializeRandomGenerator(1977+MyID);

    Pall.UnitLength = 1;
    Pall.UnitTime = 1; // Velocity unit is set to 1km/s
    Pall.UnitMass = 1;
    Pall.TCMB = 0.0;

    Pall.Lbox[0] = 1.e0;
    Pall.Lbox[1] = 0.5;
    Pall.Lbox[2] = 0.5;
    Pall.Lboxh[0] = 0.5*Pall.Lbox[0];
    Pall.Lboxh[1] = 0.5*Pall.Lbox[1];
    Pall.Lboxh[2] = 0.5*Pall.Lbox[2];
    Pall.BoxCenter[0] = 0.0;
    Pall.BoxCenter[1] = 0.25;
    Pall.BoxCenter[2] = 0.25;

    Pall.Gamma = 1.4;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;
#ifdef USE_VARIABLE_ALPHA //{
    SetVariableViscosityParameters(0.1,1.0,1.0,0.1);
#else // USE_VARIABLE_ALPHA //}//{
    SetConstantViscosityParameters(1.0);
#endif // USE_VARIABLE_ALPHA //}

    Pall.GravConst = GetUnitGravitationalConstant();
    Pall.DegreeOfFreedom = 3.0/2.0;
    Pall.HeliumWeight = 0.24;
    Pall.MeanMolecularWeight = 0.59;
    Pall.FrozenRedshift = 9.e0;
    // Pall.ConvertTtoU = GetUnitConversionFactorTemperatureToInnerEnergy();
    // Pall.ConvertUtoT = GetUnitConversionFactorInnerEnergyToTemperature();
    Pall.ConvertTtoU = 1.e-3; 
    Pall.ConvertUtoT = 1.e+3; 
    Pall.ConvertDensityToCGS = GetUnitDensityCGS();
    Pall.ConvertNumberDensityToCGS = GetUnitNumberDensityCGS();

    int NLeft = Number;
    int NRight = Number;
    int Ntotal = CUBE(NLeft)+CUBE(NRight);

    // Allocate Particle Data Structures.
    int AllocationSize = 0;
    for(int i=0;i<Ntotal;i++)
        if(i%NProcs == MyID)
            AllocationSize ++;

    dprintlmpi(AllocationSize);

    GenerateStructPbody(AllocationSize);
    GenerateStructPhydro(AllocationSize);
    for(int i=0;i<AllocationSize;i++){
        Pbody[i]->Baryon = (void *)(Phydro[i]);
        Phydro[i]->Body = Pbody[i];
    }

    double dx = 0.5/((double)Number);

    int mycount = 0;

    double mass = 0.5*0.5*1.0*2.4e4/((double)Ntotal);
    if(mode == 0){
    for(int l=0;l<2;l++){
    for(int i=0;i<NLeft;i++){
        for(int j=0;j<NLeft;j++){
            for(int k=0;k<NLeft;k++){
                PhydroBody(mycount)->Active = ON;
                PhydroBody(mycount)->Use = ON;
                PhydroBody(mycount)->Type = TypeHydro;
                PhydroBody(mycount)->GlobalID = mycount;

                if(l == 0){
                    PhydroPos(mycount)[0] = (i+0.5)*dx-0.5;
                } else {
                    PhydroPos(mycount)[0] = (i+0.5)*dx;
                }
                PhydroPos(mycount)[1] = (j+0.5)*dx;
                PhydroPos(mycount)[2] = (k+0.5)*dx;

                Phydro[mycount]->Kernel = 2.0*dx;
                Phydro[mycount]->Rho = 2.4e+4;
                Phydro[mycount]->U = 50.0;


                double t = 0.000001;
                double Kappa = 1.0/(Phydro[mycount]->Rho*Pall.ConvertTtoU);
                double X = PhydroPos(mycount)[0]/2.0/sqrt(Kappa*t);
                double U = (50.0/2.0)*erfc(X)+(50.0/T_JUMP/2.0)*erfc(-X);
                Phydro[mycount]->U = U;
                PhydroPos(mycount)[0] += 0.5;


#ifdef EVALUATE_KERNEL_BY_ITERATION
                Phydro[mycount]->Kernel = 1.e0/((double)Number);
#endif
                PhydroBody(mycount)->Mass = PhydroMass(mycount) = mass;
                //PhydroPos(mycount)[1] = PhydroPos(mycount)[2] = 0.e0;
                PhydroVel(mycount)[0] = PhydroVel(mycount)[1] = PhydroVel(mycount)[2] = 0.e0;
                Phydro[mycount]->HydroAcc[0] = Phydro[mycount]->HydroAcc[1] = Phydro[mycount]->HydroAcc[2] = 0.e0;
                Phydro[mycount]->Use = ON; 
                PhydroBody(mycount)->Type = TypeHydro;
                PhydroBody(mycount)->Eps = 1.0;
                mycount ++;

            }
        }
    }
    }
    } else {
        for(int i=0;i<Ntotal;i++){
            PhydroBody(mycount)->Active = ON;
            PhydroBody(mycount)->Use = ON;
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->GlobalID = mycount;

            PhydroPos(mycount)[0] = gsl_rng_uniform(RandomGenerator)-0.5;
            PhydroPos(mycount)[1] = gsl_rng_uniform(RandomGenerator);
            PhydroPos(mycount)[2] = gsl_rng_uniform(RandomGenerator);

            Phydro[mycount]->Kernel = 2.0*dx;
            Phydro[mycount]->Rho = 2.4e+4;

            double t = 0.000001;
            double Kappa = 1.0/(Phydro[mycount]->Rho*Pall.ConvertTtoU);
            double X = PhydroPos(mycount)[0]/2.0/sqrt(Kappa*t);
            double U = (50.0/2.0)*erfc(X)+(50.0/T_JUMP/2.0)*erfc(-X);
            Phydro[mycount]->U = U;
            PhydroPos(mycount)[0] += 0.5;


#ifdef EVALUATE_KERNEL_BY_ITERATION
            Phydro[mycount]->Kernel = 1.e0/((double)Number);
#endif
            PhydroBody(mycount)->Mass = PhydroMass(mycount) = 0.5*0.5*1.0*2.4e4/((double)Ntotal);
            //PhydroPos(mycount)[1] = PhydroPos(mycount)[2] = 0.e0;
            PhydroVel(mycount)[0] = PhydroVel(mycount)[1] = PhydroVel(mycount)[2] = 0.e0;
            Phydro[mycount]->HydroAcc[0] = Phydro[mycount]->HydroAcc[1] = Phydro[mycount]->HydroAcc[2] = 0.e0;
            Phydro[mycount]->Use = ON; 
            PhydroBody(mycount)->Type = TypeHydro;
            PhydroBody(mycount)->Eps = 1.0;
            mycount ++;
        }
    }

    Pall.Ntotal = Pall.Nhydro = mycount;
    Pall.Ntotal_t = Pall.Nhydro_t = Ntotal;

    ActivateAllparticles();

    FILE *fp;
    char fcloud[MaxCharactersInLine];
    Snprintf(fcloud,"_hydro.%03d",MyID);
    FileOpen(fp,fcloud,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%d %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroBody(i)->Pos[0],PhydroBody(i)->Pos[1],PhydroBody(i)->Pos[2],
                Phydro[i]->Rho,Phydro[i]->Kernel,Phydro[i]->U);
    }
    fclose(fp);
    fflush(NULL);
    if(MPIGetMyID() == MPI_ROOT_RANK){
        system("cat _hydro.* > hydro.dat");
        fflush(NULL);
        system("rm -rf ./_hydro.*");
        fflush(NULL);
    }
    fflush(NULL);


    Pall.RunStatus = NewSimulation;
    Pall.AdaptiveSofteningFactor = 1.0;

    Pall.Ns = 32;
    Pall.Npm = 2;
    // Pall.Ns = 128;
    // Pall.Npm = 8;

    Pall.TEnd = 0.024;
    Pall.TCurrent = 0.e0;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    Pall.OutPutFileNumber = 0;
    Pall.OutPutInterval = 0.1*Pall.TEnd;
    MakeDir("./data");

    strcpy(Pall.ASCIIFileName,"./data/TC3D.ASCII");
    strcpy(Pall.BaseFileName,"./data/TC3D");
    strcpy(Pall.RestartFileName,"./data/TC3D.dump");


    return;
}
#endif // TASK_TEST_3D_THERMAL_CONDUCTIVITY //}
