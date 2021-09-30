#include "config.h"
#include "StructureOperation.h"
#include "ParallelOperation.h"
#include "Logs.h"

#define E_BARYON(x)    ( pow((x)/(1.0e-5)*CUBE(5.0e-4),1.0/3.0) )
#define E_CDM(x)    ( pow((x)/(1.0e-5)*CUBE(5.0e-4),1.0/3.0) )

void ReadCosmologicalData(void){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));

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


    for(int i=0;i<NProcs;i++){
        if(MyID == i){
            int number,number1,number2;
            double Pos[NProcs][3],Vel[NProcs][3],Mass[NProcs];
            FILE *fp1,*fp2;

            FileOpen(fp1,"./tophat.init.bin.baryon","rb");
            FileOpen(fp2,"./tophat.init.bin.cdm","rb");

            /** Get Inital Parameter **/
            fread(&number1,sizeof(int),1,fp1);
            fprintf(stderr,"Current Number of SPH Particles = %d\n",number1);
            fseek(fp1,sizeof(int),SEEK_CUR);
            fseek(fp1,sizeof(double),SEEK_CUR);
            fread(&(Pall.InitialRedshift),sizeof(double),1,fp1);
            fseek(fp1,sizeof(double),SEEK_CUR);

            /*** Get Cosmological Parameter ***/
            fread(&(Pall.hubble),sizeof(double),1,fp1);
            fread(&(Pall.OmegaM),sizeof(double),1,fp1);
            fread(&(Pall.OmegaL),sizeof(double),1,fp1);
            fread(&(Pall.OmegaB),sizeof(double),1,fp1);
            Pall.hubble *= 0.01;
            Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

            fread(&number2,sizeof(int),1,fp2);
            fseek(fp2,sizeof(int),SEEK_CUR);
            fseek(fp2,sizeof(double)*7,SEEK_CUR);

            number = number1+number2;

            fprintf(stderr,"Cosmological Paramters\n");
            fprintf(stderr,"(Hubble,OmegaM,OmegaL,RedShift) = (%g,%g,%g,%g)\n",
                Pall.Hubble,Pall.OmegaM,Pall.OmegaL,Pall.InitialRedshift);
            fprintf(stderr,"Current Number of SPH Particles = %d\n",number1);
            fprintf(stderr,"Current Number of DM Particles = %d\n",number2);

            int counter1 = 0;  //Number of baryon particles
            for(int k=0;k<number1;k+=NProcs){
                int Element = MIN(number1-k,NProcs);
                if(MyID < Element)
                    counter1 ++;
            }
            int counter2 = 0; // Number of cold dark matter particles
            for(int k=0;k<number2;k+=NProcs){
                int Element = MIN(number2-k,NProcs);
                if(MyID < Element)
                    counter2 ++;
            }
            fprintf(stderr,"counter1 = %d, counter2 = %d\n",counter1,counter2);

            GenerateStructPbody(counter1+counter2);
            GenerateStructPhydro(counter1);
            fprintf(stderr,"PbodySize = %d, PhydroSize = %d\n",PbodySize,PhydroSize);

            // Connect Body and Hydro locally.
            for(int k=0;k<counter1;k++){
                Pbody[k] = PbodyElements+k;
                Phydro[k] = PhydroElements+k;
                Pbody[k]->Baryon = (void *)(Phydro[k]);
                Phydro[k]->Body = Pbody[k];
            }


            int counter = 0;
            for(int k=0;k<number1;k+=NProcs){
                int Element = MIN(number1-k,NProcs);
                fread(Pos,sizeof(double),3*Element,fp1);
                if(MyID < Element){
                    PbodyElements[counter].Pos[0] = Pos[MyID][0];
                    PbodyElements[counter].Pos[1] = Pos[MyID][1];
                    PbodyElements[counter].Pos[2] = Pos[MyID][2];
                    PbodyElements[counter].GlobalID = k+MyID;
                    counter++;
                }
            }
            for(int k=0;k<number2;k+=NProcs){
                int Element = MIN(number2-k,NProcs);
                fread(Pos,sizeof(double),3*Element,fp2);
                if(MyID < Element){
                    PbodyElements[counter].Pos[0] = Pos[MyID][0];
                    PbodyElements[counter].Pos[1] = Pos[MyID][1];
                    PbodyElements[counter].Pos[2] = Pos[MyID][2];
                    PbodyElements[counter].GlobalID = number1+k+MyID;
                    counter++;
                }
            }

            counter = 0;
            for(int k=0;k<number1;k+=NProcs){
                int Element = MIN(number1-k,NProcs);
                fread(Vel,sizeof(double),3*Element,fp1);
                if(MyID < Element){
                    PbodyElements[counter].Vel[0] = Vel[MyID][0];
                    PbodyElements[counter].Vel[1] = Vel[MyID][1];
                    PbodyElements[counter].Vel[2] = Vel[MyID][2];
                    counter++;
                }
            }
            for(int k=0;k<number2;k+=NProcs){
                int Element = MIN(number2-k,NProcs);
                fread(Vel,sizeof(double),3*Element,fp2);
                if(MyID < Element){
                    PbodyElements[counter].Vel[0] = Vel[MyID][0];
                    PbodyElements[counter].Vel[1] = Vel[MyID][1];
                    PbodyElements[counter].Vel[2] = Vel[MyID][2];
                    counter++;
                }
            }
            
            counter = 0;
            for(int k=0;k<number1;k+=NProcs){
                int Element = MIN(number1-k,NProcs);
                fread(Mass,sizeof(double),Element,fp1);
                if(MyID < Element){
                    PbodyElements[counter].Mass = Mass[MyID];
                    counter++;
                }
            }
            for(int k=0;k<number2;k+=NProcs){
                int Element = MIN(number2-k,NProcs);
                fread(Mass,sizeof(double),Element,fp2);
                if(MyID < Element){
                    PbodyElements[counter].Mass = Mass[MyID];
                    counter++;
                }
            }
    
            fclose(fp1);
            fclose(fp2);

            Pall.Ntotal = counter1+counter2;
            Pall.Nhydro = counter1;
            Pall.NDM = counter2;

            for(int k=0;k<Pall.Nhydro;k++){
                PhydroElements[k].Use = ON;
                PhydroElements[k].Active = ON;
            }
            for(int k=0;k<Pall.Ntotal;k++){
                PbodyElements[k].Use = ON;
                PbodyElements[k].Active = ON;
            }

            ReConnectPointers();
            UpdateTotalNumber();

            //double EPS_BARYON = E_BARYON(PhydroMass(0));
            double EPS_BARYON = E_BARYON(PbodyElements[0].Mass);
            for(int k=0;k<Pall.Nhydro;k++){
                PhydroActive(k) = ON;
                Phydro[k]->U = Pall.ConvertTtoU*Pall.TCMB*(1.e0+Pall.InitialRedshift);
                PhydroBody(k)->Eps = EPS_BARYON;
                Phydro[k]->Kernel = PhydroBody(k)->Eps;
                PhydroBody(k)->Type = TypeHydro;

#if (UseSFModelSpawn) 
                PbodyHydro(k)->SpawnMass = Pbody[k]->Mass/(double)MaxSpawnTimes;
#endif
                PbodyHydro(k)->ZII = 0.e0;
                PbodyHydro(k)->ZIa = 0.e0;
            }

            //double EPS_CDM = E_CDM(Pbody[counter1]->Mass);
            double EPS_CDM = EPS_BARYON;
            for(int k=counter1;k<Pall.Ntotal;k++){
                Pbody[k]->Active = ON;
                Pbody[k]->Eps = EPS_CDM;
                Pbody[k]->Type = TypeDM;
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


    // print out softening length
    fprintf(stderr,"Baryonic softening length = %g kpc.\n",1000*E_BARYON(PhydroMass(0)));

#if 0
    // Plot Particle Parameters
    {
    FILE *fp;
    char Fname[MaxCharactersInLine];
    sprintf(Fname,"Cosmological.Init.Hydro.%02d.%02d.data",MPIGetMyID(),MPIGetNumProcs());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2]);
    }
    fclose(fp);

    sprintf(Fname,"Cosmological.Init.DM.%02d.%02d.data",MPIGetMyID(),MPIGetNumProcs());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            fprintf(fp,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
        }
    }
    fclose(fp);
    }
#endif

    Pall.RunStatus = NewSimulation;

    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = CalcZtoT(0.e0);
    Pall.TCurrent = CalcZtoT(Pall.InitialRedshift);
    Pall.Redshift = Pall.InitialRedshift;
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
            Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.e0;
    Pall.ViscousS = 100.e0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    UpdateAdaptiveSofteningFactor();

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    InitializeRandomGenerator(1977+MPIGetMyID());

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = 0.01*Pall.UnitTime/GIGAYEAR_CGS; 
    //Pall.OutPutInterval = 0.01*GIGAYEAR_CGS/Pall.UnitTime; 
    //Pall.OutPutInterval = GIGAYEAR_CGS/Pall.UnitTime; 
    Pall.OutPutInterval = 0.1*GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/Cosmological.ASCII");
    strcpy(Pall.BaseFileName,"./data/Cosmological");
    strcpy(Pall.RestartFileName,"./data/Cosmological.dump");
    InitLogFiles();

    return;
}

void ReadCosmologicalHalfwayData(char *fname){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));

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

    InitializeRandomGenerator(1977+MPIGetMyID());

    struct StructHydroBasket{
        double Pos[3];
        double Vel[3];
        double Mass;
        double Rho;
        double Kernel;
        double U;
        unsigned long int GlobalID;
    } HydroBasket;

    struct StructStarBasket{
        double Pos[3];
        double Vel[3];
        double Mass;
        double Age;
        unsigned long int GlobalID;
    } StarBasket;

    struct StructDMBasket{
        double Pos[3];
        double Vel[3];
        double Mass;
        double Kernel;
        unsigned long int GlobalID;
    } DMBasket;


    Pall.hubble = 50;
    Pall.OmegaM = 1;
    Pall.OmegaL = 0;
    Pall.OmegaB = 0.1;
    Pall.hubble *= 0.01;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

    fprintf(stderr,"Cosmological Paramters\n");
    fprintf(stderr,"(Hubble,OmegaM,OmegaL) = (%g,%g,%g)\n",
        Pall.Hubble,Pall.OmegaM,Pall.OmegaL);

    FILE *fp;
    FileOpen(fp,fname,"r");
    fread(&Pall.Nhydro_t,sizeof(int),1,fp);
    fread(&Pall.Nstars_t,sizeof(int),1,fp);
    fread(&Pall.NDM_t,sizeof(int),1,fp);
    Pall.Ntotal_t = Pall.Nhydro_t + Pall.Nstars_t + Pall.NDM_t;

    fread(&Pall.TCurrent,sizeof(double),1,fp);
    fread(&Pall.TStepTotal,sizeof(int),1,fp);

    int count_hydro = 0;
    for(int i=0;i<Pall.Nhydro_t;i++)
        if(i%NProcs == MyID)
            count_hydro ++;

    int count_star = 0;
    for(int i=0;i<Pall.Nstars_t;i++)
        if(i%NProcs == MyID)
            count_star ++;

    int count_dm = 0;
    for(int i=0;i<Pall.NDM_t;i++)
        if(i%NProcs == MyID)
            count_dm ++;

    int count_hydro_t,count_star_t,count_dm_t;
    MPI_Allreduce(&count_hydro,&count_hydro_t,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&count_star,&count_star_t,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&count_dm,&count_dm_t,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    assert(Pall.Nhydro_t == count_hydro_t);
    assert(Pall.Nstars_t == count_star_t);
    assert(Pall.NDM_t == count_dm_t);

    Pall.Nhydro = count_hydro;
    Pall.Nstars = count_star;
    Pall.NDM = count_dm;
    Pall.Ntotal = count_hydro + count_star + count_dm;

    fprintf(stderr,"Number of Total Particles : %ld %ld %ld\n",Pall.Nhydro_t,Pall.Nstars_t,Pall.NDM_t);
    fprintf(stderr,"Number of Particles : %ld %ld %ld\n",Pall.Nhydro,Pall.Nstars,Pall.NDM);


    GenerateStructPbody(Pall.Ntotal);
    GenerateStructPhydro(Pall.Nhydro);
    GenerateStructPstar(Pall.Nstars);

    // Phydro
    count_hydro = 0;
    for(int i=0;i<Pall.Nhydro_t;i++){
        fread(&HydroBasket,sizeof(struct StructHydroBasket),1,fp);
        if(i%NProcs == MyID){
            PbodyElements[count_hydro].Pos[0] = HydroBasket.Pos[0];
            PbodyElements[count_hydro].Pos[1] = HydroBasket.Pos[1];
            PbodyElements[count_hydro].Pos[2] = HydroBasket.Pos[2];

            PbodyElements[count_hydro].Vel[0] = HydroBasket.Vel[0];
            PbodyElements[count_hydro].Vel[1] = HydroBasket.Vel[1];
            PbodyElements[count_hydro].Vel[2] = HydroBasket.Vel[2];

            PbodyElements[count_hydro].Mass = HydroBasket.Mass;

            PbodyElements[count_hydro].GlobalID = HydroBasket.GlobalID;

            PhydroElements[count_hydro].Rho = HydroBasket.Rho;
            PhydroElements[count_hydro].Kernel = HydroBasket.Kernel;
            PhydroElements[count_hydro].U = HydroBasket.U;

            PbodyElements[count_hydro].Baryon = PhydroElements+count_hydro;
            PhydroElements[count_hydro].Body = PbodyElements+count_hydro;

            // others for run
            PbodyElements[count_hydro].Type = TypeHydro;
            PbodyElements[count_hydro].Active = ON;

#if (UseSFModelSpawn) 
            PhydroElements[count_hydro].SpawnMass = PbodyElements[count_hydro].Mass/(double)MaxSpawnTimes;
#endif
            PhydroElements[count_hydro].ZII = 0.e0;
            PhydroElements[count_hydro].ZIa = 0.e0;

            PhydroElements[count_hydro].Use = ON;

            count_hydro ++;
        }
    }

    // Pstar
    count_star = 0;
    for(int i=0;i<Pall.Nstars_t;i++){
        fread(&StarBasket,sizeof(struct StructStarBasket),1,fp);
        if(i%NProcs == MyID){
            PbodyElements[Pall.Nhydro+count_star].Pos[0] = StarBasket.Pos[0];
            PbodyElements[Pall.Nhydro+count_star].Pos[1] = StarBasket.Pos[1];
            PbodyElements[Pall.Nhydro+count_star].Pos[2] = StarBasket.Pos[2];

            PbodyElements[Pall.Nhydro+count_star].Vel[0] = StarBasket.Vel[0];
            PbodyElements[Pall.Nhydro+count_star].Vel[1] = StarBasket.Vel[1];
            PbodyElements[Pall.Nhydro+count_star].Vel[2] = StarBasket.Vel[2];

            PbodyElements[Pall.Nhydro+count_star].Mass = StarBasket.Mass;
            PbodyElements[Pall.Nhydro+count_star].GlobalID = StarBasket.GlobalID;

            PstarElements[count_star].FormationTime = StarBasket.Age;

            PbodyElements[Pall.Nhydro+count_star].Baryon = PstarElements+count_star;
            PstarElements[count_star].Body = PbodyElements+Pall.Nhydro+count_star;

            // others for run
            PbodyElements[Pall.Nhydro+count_star].Type = TypeStar;
            PbodyElements[Pall.Nhydro+count_star].Active = ON;

#define TSNIIEND 3.e+7 // 30 Myr average.
            if(Pall.TCurrent - PstarElements[count_star].FormationTime > TSNIIEND*YEAR_CGS/Pall.UnitTime){
                PstarElements[count_star].TypeII = true;
            } else {
                double p = (Pall.TCurrent - PstarElements[count_star].FormationTime)/(TSNIIEND*YEAR_CGS/Pall.UnitTime);
                if(gsl_rng_uniform(RandomGenerator)<p){ 
                    PstarElements[count_star].TypeII = true;
                } else {
                    PstarElements[count_star].TypeII = false;
                }
            }
            PstarElements[count_star].TypeIa = false;

            PstarElements[count_star].ParentGlobalID = PbodyElements[Pall.Nhydro+count_star].GlobalID;
            PstarElements[count_star].NthChildren = 1;
            PstarElements[count_star].Mass = PbodyElements[Pall.Nhydro+count_star].Mass;
            PstarElements[count_star].InitialMass = PbodyElements[Pall.Nhydro+count_star].Mass;
            PstarElements[count_star].Z = 0.e0;
            PstarElements[count_star].ZII = 0.e0;
            PstarElements[count_star].ZIa = 0.e0;

            PstarElements[count_star].Use = ON;
            count_star ++;
        }
    }

    count_dm = 0;
    for(int i=0;i<Pall.NDM_t;i++){
        fread(&DMBasket,sizeof(struct StructDMBasket),1,fp);
        if(i%NProcs == MyID){
            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Pos[0] = DMBasket.Pos[0];
            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Pos[1] = DMBasket.Pos[1];
            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Pos[2] = DMBasket.Pos[2];

            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Vel[0] = DMBasket.Vel[0];
            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Vel[1] = DMBasket.Vel[1];
            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Vel[2] = DMBasket.Vel[2];

            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Mass = DMBasket.Mass;
            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].GlobalID = DMBasket.GlobalID;

            // others for run
            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Type = TypeDM;
            PbodyElements[Pall.Nhydro+Pall.Nstars+count_dm].Active = ON;

            count_dm ++;
        }
    }

    fclose(fp);

    // use on set softening

    double eps = E_BARYON(((StructPbodyptr)(PhydroElements[0].Body))->Mass);
    for(int i=0;i<Pall.Ntotal;i++){
        PbodyElements[i].Use = ON;
        PbodyElements[i].Eps = eps;
    }
    ReConnectPointers();
    UpdateTotalNumber();
    Pall.NActives_t = Pall.Ntotal_t;

    fprintf(stderr,"Number of Total Particles : %ld %ld %ld\n",Pall.Nhydro_t,Pall.Nstars_t,Pall.NDM_t);
    fprintf(stderr,"Number of Particles : %ld %ld %ld\n",Pall.Nhydro,Pall.Nstars,Pall.NDM);

    // print out softening length
    fprintf(stderr,"Baryonic softening length = %g kpc.\n",1000*E_BARYON(PhydroMass(0)));

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }


#if 0
    // Plot Particle Parameters
    {
    FILE *fp;
    char Fname[MaxCharactersInLine];
    sprintf(Fname,"Cosmological.Init.Hydro.%02d.%02d.data",MPIGetMyID(),MPIGetNumProcs());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2]);
    }
    fclose(fp);

    sprintf(Fname,"Cosmological.Init.DM.%02d.%02d.data",MPIGetMyID(),MPIGetNumProcs());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            fprintf(fp,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
        }
    }
    fclose(fp);
    }

    exit(0);
#endif

    Pall.RunStatus = NewSimulation;

    Pall.Ns = 32;
    Pall.Npm = 2;

    Pall.TEnd = CalcZtoT(0.e0);
    //Pall.TCurrent = CalcZtoT(Pall.InitialRedshift);
    //Pall.Redshift = Pall.InitialRedshift;
    Pall.Redshift = CalcTtoZ(Pall.TCurrent);
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
            Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

#ifdef USE_VARIABLE_ALPHA
    Pall.ViscousAlphaMin = 0.01;
    Pall.ViscousAlphaMax = 1.e0;
    Pall.ViscousS = 100.e0;
    Pall.ViscousL = 5.0;
#endif // USE_VARIABLE_ALPHA

    UpdateAdaptiveSofteningFactor();

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    Pall.OutPutFileNumber = 0;
    //Pall.OutPutInterval = 0.01*Pall.UnitTime/GIGAYEAR_CGS; 
    Pall.OutPutInterval = 0.01*GIGAYEAR_CGS/Pall.UnitTime; 

    MakeDir("./data");
    strcpy(Pall.ASCIIFileName,"./data/Cosmological.ASCII");
    strcpy(Pall.BaseFileName,"./data/Cosmological");
    strcpy(Pall.RestartFileName,"./data/Cosmological.dump");

    return;
}

void ReadOldCosmologicalData(const int FileID){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));

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


    for(int i=0;i<NProcs;i++){
        if(MyID == i){
            int skip;
            int number,number_hydro,number_dm,int_dummy;
            double double_dummy;
            double Pos[NProcs][3],Vel[NProcs][3],Mass[NProcs];
            FILE *fp;
            char fname[MaxCharactersInLine];
            sprintf(fname,"./out%03d",FileID);
            FileOpen(fp,fname,"rb");

            /** Get Inital Parameter **/
            fread(&skip,sizeof(int),1,fp); // skip
            fread(&number,sizeof(int),1,fp); // Ntotal_t
            fprintf(stderr,"Current Number of SPH Particles = %d\n",number); 
            fread(&number_hydro,sizeof(int),1,fp); // Nhydro_t
            fread(&int_dummy,sizeof(int),1,fp);    // dummy
            fread(&int_dummy,sizeof(int),1,fp);    // dummy
            fread(&number_dm,sizeof(int),1,fp);    // Ndm_t
            fread(&(Pall.hubble),sizeof(double),1,fp);  // hubble 
            fread(&(Pall.OmegaM),sizeof(double),1,fp);  // omega M
            fread(&(Pall.OmegaL),sizeof(double),1,fp);  // omega L

            fread(&(Pall.TCurrent),sizeof(double),1,fp);  // t
            fread(&double_dummy,sizeof(double),1,fp);  // tnext
            fread(&int_dummy,sizeof(int),1,fp);    // step

            number = number_hydro+number_dm;

            fprintf(stderr,"Cosmological Paramters\n");
            fprintf(stderr,"(Hubble,OmegaM,OmegaL,RedShift) = (%g,%g,%g,%g)\n",
                Pall.Hubble,Pall.OmegaM,Pall.OmegaL,Pall.InitialRedshift);
            fprintf(stderr,"Current Number of SPH Particles = %d\n",number_hydro);
            fprintf(stderr,"Current Number of DM Particles = %d\n",number_dm);

            int counter = 0; // Number of cold dark matter particles
            for(int k=0;k<number;k+=NProcs){
                int Element = MIN(number-k,NProcs);
                if(MyID < Element)
                    counter ++;
            }

            GenerateStructPbody(counter);
            fprintf(stderr,"PbodySize = %d\n",PbodySize);

            // Connect Body and Hydro locally.
            for(int k=0;k<counter;k++)
                Pbody[k] = PbodyElements+k;

            fseek(fp,number*sizeof(int),SEEK_CUR); // skip status

            counter = 0;
            for(int k=0;k<number;k+=NProcs){
                int Element = MIN(number-k,NProcs);
                fread(Pos,sizeof(double),3*Element,fp);
                if(MyID < Element){
                    PbodyElements[counter].Pos[0] = Pos[MyID][0];
                    PbodyElements[counter].Pos[1] = Pos[MyID][1];
                    PbodyElements[counter].Pos[2] = Pos[MyID][2];
                    PbodyElements[counter].GlobalID = k+MyID;
                    counter++;
                }
            }

            counter = 0;
            for(int k=0;k<number;k+=NProcs){
                int Element = MIN(number-k,NProcs);
                fread(Vel,sizeof(double),3*Element,fp);
                if(MyID < Element){
                    PbodyElements[counter].Vel[0] = Vel[MyID][0];
                    PbodyElements[counter].Vel[1] = Vel[MyID][1];
                    PbodyElements[counter].Vel[2] = Vel[MyID][2];
                    counter++;
                }
            }
            
            counter = 0;
            for(int k=0;k<number;k+=NProcs){
                int Element = MIN(number-k,NProcs);
                fread(Mass,sizeof(double),Element,fp);
                if(MyID < Element){
                    PbodyElements[counter].Mass = Mass[MyID];
                    counter++;
                }
            }
            fclose(fp);
    
            Pall.Ntotal = counter;
            Pall.Nhydro = 0;
            Pall.NDM = counter;

            for(int k=0;k<Pall.Ntotal;k++)
                PbodyElements[k].Use = ON;

            ReConnectPointers();
            UpdateTotalNumber();

            double MassMin = Pbody[0]->Mass;
            for(int k=1;k<Pall.Ntotal;k++)
                MassMin = fmin(MassMin,Pbody[k]->Mass);

            double EPS_CDM = E_BARYON(MassMin);
            for(int k=0;k<Pall.Ntotal;k++){
                Pbody[k]->Active = ON;
                Pbody[k]->Eps = EPS_CDM;
                Pbody[k]->Type = TypeDM;
                Pbody[k]->InteractionList = 1;
            }
            fprintf(stderr,"Baryonic softening length = %g kpc.\n",1000*E_BARYON(MassMin));
        }
    }

#if 0
    // Plot Particle Parameters
    {
    FILE *fp;
    char Fname[MaxCharactersInLine];
    sprintf(Fname,"Cosmological.Init.Hydro.%02d.%02d.data",MPIGetMyID(),MPIGetNumProcs());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
                PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
                PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2]);
    }
    fclose(fp);

    sprintf(Fname,"Cosmological.Init.DM.%02d.%02d.data",MPIGetMyID(),MPIGetNumProcs());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            fprintf(fp,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
        }
    }
    fclose(fp);
    }
#endif

    /*
    {
    FILE *fp;
    char Fname[MaxCharactersInLine];
    sprintf(Fname,"Cosmological.Init.%02d.%02d.data",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
    }
    fclose(fp);
    }
    */

    Pall.RunStatus = NewSimulation;

    //Pall.Ns = 50;
    //Pall.Npm = 5;

    Pall.hubble = 0.01*50;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

    Pall.TEnd = CalcZtoT(0.e0);
    //Pall.TCurrent = CalcZtoT(Pall.InitialRedshift);
    Pall.Redshift = CalcTtoZ(Pall.TCurrent);
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
            Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0*Pall.HydroAlpha;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    UpdateAdaptiveSofteningFactor();

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
            Pall.AdaptiveSofteningFactor);

    InitLogFiles();

    //MPI_Finalize();
    //exit(EXIT_SUCCESS);

    return;
}

void ReadOldCosmologicalDataFull(const int FileID){

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();

    memset(&Pall,0,sizeof(struct StructPall));

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

    for(int i=0;i<NProcs;i++){
        if(MyID == i){
            int skip;
            int number,number_hydro,number_dm,int_dummy;
            double double_dummy;
            double Pos[NProcs][3],Vel[NProcs][3],Mass[NProcs];
            FILE *fp;
            char fname[MaxCharactersInLine];
            sprintf(fname,"./out%03d",FileID);
            FileOpen(fp,fname,"rb");

            /** Get Inital Parameter **/
            fread(&skip,sizeof(int),1,fp); // skip
            fread(&number,sizeof(int),1,fp); // Ntotal_t
            fprintf(stderr,"Current Number of SPH Particles = %d\n",number);
            fread(&number_hydro,sizeof(int),1,fp); // Nhydro_t
            fread(&int_dummy,sizeof(int),1,fp);    // dummy
            fread(&int_dummy,sizeof(int),1,fp);    // dummy
            fread(&number_dm,sizeof(int),1,fp);    // Ndm_t
            fread(&(Pall.hubble),sizeof(double),1,fp);  // hubble
            fread(&(Pall.OmegaM),sizeof(double),1,fp);  // omega M
            fread(&(Pall.OmegaL),sizeof(double),1,fp);  // omega L

            fread(&(Pall.TCurrent),sizeof(double),1,fp);  // t
            fread(&double_dummy,sizeof(double),1,fp);  // tnext
            fread(&int_dummy,sizeof(int),1,fp);    // step

            number = number_hydro+number_dm;

            fprintf(stderr,"Cosmological Paramters\n");
            fprintf(stderr,"(Hubble,OmegaM,OmegaL,RedShift) = (%g,%g,%g,%g)\n",
            Pall.Hubble,Pall.OmegaM,Pall.OmegaL,Pall.InitialRedshift);
            fprintf(stderr,"Current Number of SPH Particles = %d\n",number_hydro);
            fprintf(stderr,"Current Number of DM Particles = %d\n",number_dm);


            int counter = 0; // Number of cold dark matter particles
            for(int k=0;k<number;k+=NProcs){
                int Element = MIN(number-k,NProcs);
                if(MyID < Element)
                    counter ++;
            }
            GenerateStructPbody(counter);
            fprintf(stderr,"PbodySize = %d\n",PbodySize);

            int counter1 = 0; // Number of cold dark matter particles
            for(int k=0;k<number_hydro;k+=NProcs){
                int Element = MIN(number_hydro-k,NProcs);
                if(MyID < Element)
                    counter1 ++;
            }
            GenerateStructPhydro(counter1);
            fprintf(stderr,"PhydroSize = %d\n",PhydroSize);

            for(int k=0;k<counter;k++)
                Pbody[k] = PbodyElements+k;

            // Connect Body and Hydro locally.
            for(int k=0;k<counter1;k++){
                Pbody[k] = PbodyElements+k;
                Phydro[k] = PhydroElements+k;
                Pbody[k]->Baryon = (void *)(Phydro[k]);
                Phydro[k]->Body = Pbody[k];
            }

            fseek(fp,number*sizeof(int),SEEK_CUR); // skip status


            counter = 0;
            for(int k=0;k<number;k+=NProcs){
                int Element = MIN(number-k,NProcs);
                fread(Pos,sizeof(double),3*Element,fp);
                if(MyID < Element){
                    PbodyElements[counter].Pos[0] = Pos[MyID][0];
                    PbodyElements[counter].Pos[1] = Pos[MyID][1];
                    PbodyElements[counter].Pos[2] = Pos[MyID][2];
                    PbodyElements[counter].GlobalID = k+MyID;
                    counter++;
                }
            }

            counter = 0;
            for(int k=0;k<number;k+=NProcs){
                int Element = MIN(number-k,NProcs);
                fread(Vel,sizeof(double),3*Element,fp);
                if(MyID < Element){
                    PbodyElements[counter].Vel[0] = Vel[MyID][0];
                    PbodyElements[counter].Vel[1] = Vel[MyID][1];
                    PbodyElements[counter].Vel[2] = Vel[MyID][2];
                    counter++;
                }
            }

            counter = 0;
            for(int k=0;k<number;k+=NProcs){
                int Element = MIN(number-k,NProcs);
                fread(Mass,sizeof(double),Element,fp);
                if(MyID < Element){
                    PbodyElements[counter].Mass = Mass[MyID];
                    counter++;
                }
            }
            fclose(fp);

            Pall.Ntotal = counter;
            Pall.Nhydro = counter1;
            Pall.NDM = counter-counter1;

            for(int k=0;k<Pall.Nhydro;k++)
                PhydroElements[k].Use = ON;
            
            for(int k=0;k<Pall.Ntotal;k++)
                PbodyElements[k].Use = ON;


            for(int k=0;k<Pall.Ntotal;k++){
                if(PbodyElements[k].GlobalID < counter1){
                    PbodyElements[k].Type = TypeHydro;
                }else{
                    PbodyElements[k].Type = TypeDM;
                }
            }

            ReConnectPointers();
            UpdateTotalNumber();

            for(int k=0;k<Pall.Nhydro;k++){
                Phydro[k]->Mass = PhydroBody(k)->Mass;
            }

            double EPS_BARYON = E_BARYON(PhydroMass(0));
            for(int k=0;k<Pall.Nhydro;k++){
                PhydroActive(k) = ON;
                Phydro[k]->U = Pall.ConvertTtoU*Pall.TCMB*(1.e0+Pall.InitialRedshift);
                PhydroBody(k)->Eps = EPS_BARYON;
                Phydro[k]->Kernel = 0.1*EPS_BARYON;
                PhydroBody(k)->Type = TypeHydro;
            }
    
            gprintlmpi(EPS_BARYON);

            //double EPS_CDM = E_BARYON(MassMin);
            double EPS_CDM = EPS_BARYON;
            for(int k=0;k<Pall.Ntotal;k++){
                Pbody[k]->Active = ON;
                if(Pbody[k]->Type != TypeHydro){
                    Pbody[k]->Eps = EPS_CDM;
                    Pbody[k]->Type = TypeDM;
                }
                Pbody[k]->InteractionList = 1;
            }

            fprintf(stderr,"Baryonic softening length = %g kpc.\n",1000*E_BARYON(PhydroMass(0)));
        }
    }

    for(int i=0;i<Pall.Nhydro;i++){
        Phydro[i]->PosP[0] = PhydroBody(i)->Pos[0];
        Phydro[i]->PosP[1] = PhydroBody(i)->Pos[1];
        Phydro[i]->PosP[2] = PhydroBody(i)->Pos[2];
        Phydro[i]->Mass = PhydroBody(i)->Mass;
        Phydro[i]->Active = PhydroBody(i)->Active;
    }

#if 0
    // Plot Particle Parameters
    {
    FILE *fp;
    char Fname[MaxCharactersInLine];
    sprintf(Fname,"Cosmological.Init.Hydro.%02d.%02d.data",MPIGetMyID(),MPIGetNumProcs());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Nhydro;i++){
        fprintf(fp,"%ld %g %g %g %g %g %g\n",PhydroBody(i)->GlobalID,
            PhydroPos(i)[0],PhydroPos(i)[1],PhydroPos(i)[2],
            PhydroVel(i)[0],PhydroVel(i)[1],PhydroVel(i)[2]);
    }
    fclose(fp);

    sprintf(Fname,"Cosmological.Init.DM.%02d.%02d.data",MPIGetMyID(),MPIGetNumProcs());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
        if(Pbody[i]->Type == TypeDM){
            fprintf(fp,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
                Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
                Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
        }
    }
    fclose(fp);
    }
#endif

    /*
    {
    FILE *fp;
    char Fname[MaxCharactersInLine];
    sprintf(Fname,"Cosmological.Init.%02d.%02d.data",MPIGetNumProcs(),MPIGetMyID());
    FileOpen(fp,Fname,"w");
    for(int i=0;i<Pall.Ntotal;i++){
    fprintf(fp,"%ld %g %g %g %g %g %g\n",Pbody[i]->GlobalID,
    Pbody[i]->Pos[0],Pbody[i]->Pos[1],Pbody[i]->Pos[2],
    Pbody[i]->Vel[0],Pbody[i]->Vel[1],Pbody[i]->Vel[2]);
    }
    fclose(fp);
    }
    */

    Pall.RunStatus = NewSimulation;

    Pall.Ns = 32;
    Pall.Npm = 5;

    Pall.hubble = 0.01*50;
    Pall.Hubble = Pall.hubble*100*(1.e+5/Pall.UnitLength)/(1.0/Pall.UnitTime);

    Pall.TEnd = CalcZtoT(0.e0);
    //Pall.TCurrent = CalcZtoT(Pall.InitialRedshift);
    Pall.Redshift = CalcTtoZ(Pall.TCurrent);
    Pall.Era = Pall.EraLocal = 0.e0;
    Pall.TStepTotal = 0;

    fprintf(stderr,"Tstart = %g, Tend = %g, Zstart = %g, Zinit = %g\n",
    Pall.TCurrent,Pall.TEnd,Pall.Redshift,Pall.InitialRedshift);

    Pall.HydroAlpha = 1.0;
    Pall.HydroBeta =  2.0*Pall.HydroAlpha;
    Pall.HydroEta2 =  0.01;
    Pall.Gamma = 5.0/3.0;
    Pall.Gm1   = Pall.Gamma-1.0;
    Pall.GGm1  = Pall.Gamma*Pall.Gm1;

    UpdateAdaptiveSofteningFactor();

    fprintf(stderr,"Pall.AdaptiveSofteningFactor = %g\n",
        Pall.AdaptiveSofteningFactor);

    InitLogFiles();

    //MPI_Finalize();
    //exit(EXIT_SUCCESS);

    return;
}


