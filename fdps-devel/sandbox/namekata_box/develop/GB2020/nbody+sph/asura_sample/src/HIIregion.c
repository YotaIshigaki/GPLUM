#include "config.h"
#include "NeighborSearch.h"
#include "StellarFeedback.h"
#include "SizeDetermination.h"

#define __LOCAL__
#define __UPDATE_ONLY_ACTIVE__
#define __SKIP_HOTGAS__
//#define __EXPORT_TIMESTEP__
#define __DT_DIFF__ (MAX_K_LOCAL)
#define __PHOTON_COUNT_BASE__

#ifdef TASK_TEST_STROMGRENSPHERE
extern bool *HIIFlag; 
#endif //

#define NumberofDataSet_nLy (10)

static double LocalKernelMax = 0.e0;

/*
 * Three members are introduced into the Pstar structure.
 * 1. StromgrenRadius <- the Stromgren radius.
 * 2. Density <- the local gas density.
 * 3. HII flag <- if the flag is ``on'', the particle has its own HII region.
 *
 * One member is introduced into the Phydro structure.
 * 1. HIIFlag <- when this flag is ``on'', the value of the FUV is set to be G0_HIIREGION.
 */

extern double SNIINumber;
int NumberofLy = 0;
//double *AgeLy;
double *nLy[NumberofDataSet_nLy];
int nLyEnd[NumberofDataSet_nLy];


double AgeLy[] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,
    100,120,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,
    1800,2000,2500,3000,3500,4000,4500,5000,6000,7000,8000,9000,10000,11000,12000,13000,
    14000,15000,16000,17000,18000,19000,20000};

double nLyZ0[] = {.282E+47,.317E+47,.348E+47,.302E+47,.169E+47,.105E+47,.694E+46,.464E+46,
    .319E+46,.221E+46,.153E+46,.720E+45,.345E+45,.170E+45,.882E+44,.498E+44,.169E+44,.780E+43,
    .426E+43,.257E+43,.166E+43,.113E+43,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00};

double nLyZSolar[] = {.243E+47,.235E+47,.193E+47,.143E+47,.681E+46,.358E+46,.106E+46,.550E+45,
    .293E+45,.170E+45,.105E+45,.460E+44,.234E+44,.131E+44,.802E+43,.518E+43,.211E+43,.991E+42,
    .527E+42,.294E+42,.173E+42,.104E+42,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
    .000E+00};

static void LoadnLyData(void){

    NumberofLy = sizeof(AgeLy)/sizeof(double);
    // if(MPIGetMyID() == MPI_ROOT_RANK){
        // dprintlmpi(NumberofLy);
    // }

    nLy[0] = malloc(sizeof(double)*NumberofLy);
    for(int i=0;i<NumberofLy;i++){
        nLy[0][i] = nLyZSolar[i];
        //fprintf(stderr,"%d %g\n",i,nLy[0][i]);
    }

    for(int i=0;i<NumberofLy;i++){
        if(nLy[0][i] == 0.e0){
            nLyEnd[0] = i;
            break;
        }
    }
    //fprintf(stderr,"End %d, %g %g\n",nLyEnd[0],
            //nLy[0][nLyEnd[0]],nLy[0][nLyEnd[0]+1]);

    return ;
}


#if 0
Ejecta of massive stars: SNII (model B of Woosley & Weaver [1995]) + stellar winds                  
Initial mass function: IMF_Salpeter.dat from   0.1000 to 100.0000 solar masses                      
Evolutionary tracks: tracksZ0.0001.dat                                                              
Metallicity (mass fraction): .0001                                                                  
Evolutionary tracks: tracksZ0.0004.dat                                                              
Metallicity (mass fraction): .0004                                                                  
Evolutionary tracks: tracksZ0.004.dat                                                               
Metallicity (mass fraction): .0040                                                                  
Evolutionary tracks: tracksZ0.008.dat                                                               
Metallicity (mass fraction): .0080                                                                  
Evolutionary tracks: tracksZ0.02.dat                                                                
Metallicity (mass fraction): .0200                                                                  
Evolutionary tracks: tracksZ0.05.dat                                                                
Metallicity (mass fraction): .0500                                                                  
Evolutionary tracks: tracksZ0.1.dat                                                                 
Metallicity (mass fraction): .1000                                                                  
Fraction of close binary systems: 0.50000E-01                                                       
Initial metallicity: 0.00000E+00                                                                    
No infall                                                                                           
Type of star formation:   0                                                                         
Consistent evolution of the stellar metallicity                                                     
Mass fraction of substellar objects: 0.00000E+00                                                    
No galactic winds                                                                                   
Nebular emission                                                                                    
No extinction                                                                                       
 time nLymcont   L(Ha)    W(Ha)    L(Hb)    W(Hb)  LB/LBsol LV/LVsol   D4000
    0 .282E+47 .391E+35 .355E+04 .134E+35 .703E+03 .147E+03 .104E+03 .651E+00
    1 .317E+47 .439E+35 .359E+04 .151E+35 .713E+03 .164E+03 .117E+03 .650E+00
    2 .348E+47 .482E+35 .351E+04 .165E+35 .686E+03 .183E+03 .129E+03 .653E+00
    3 .302E+47 .418E+35 .228E+04 .143E+35 .371E+03 .215E+03 .137E+03 .711E+00
    4 .169E+47 .235E+35 .468E+03 .805E+34 .105E+03 .246E+03 .180E+03 .887E+00
    5 .105E+47 .146E+35 .118E+04 .501E+34 .167E+03 .127E+03 .716E+02 .769E+00
    6 .694E+46 .962E+34 .136E+04 .330E+34 .189E+03 .774E+02 .434E+02 .753E+00
    7 .464E+46 .644E+34 .119E+04 .221E+34 .162E+03 .583E+02 .318E+02 .763E+00
    8 .319E+46 .442E+34 .967E+03 .151E+34 .128E+03 .482E+02 .253E+02 .776E+00
    9 .221E+46 .306E+34 .778E+03 .105E+34 .101E+03 .406E+02 .207E+02 .787E+00
   10 .153E+46 .212E+34 .622E+03 .726E+33 .794E+02 .346E+02 .171E+02 .798E+00
   12 .720E+45 .999E+33 .359E+03 .343E+33 .449E+02 .273E+02 .129E+02 .819E+00
   14 .345E+45 .479E+33 .192E+03 .164E+33 .238E+02 .237E+02 .109E+02 .835E+00
   16 .170E+45 .236E+33 .107E+03 .811E+32 .133E+02 .206E+02 .942E+01 .847E+00
   18 .882E+44 .122E+33 .587E+02 .419E+32 .736E+01 .191E+02 .868E+01 .856E+00
   20 .498E+44 .690E+32 .343E+02 .237E+32 .435E+01 .181E+02 .829E+01 .864E+00
   25 .169E+44 .235E+32 .123E+02 .806E+31 .160E+01 .166E+02 .774E+01 .878E+00
   30 .780E+43 .108E+32 .612E+01 .371E+31 .808E+00 .151E+02 .712E+01 .887E+00
   35 .426E+43 .783E+31 .486E+01 .268E+31 .650E+00 .137E+02 .648E+01 .896E+00
   40 .257E+43 .591E+31 .386E+01 .203E+31 .523E+00 .128E+02 .614E+01 .903E+00
   45 .166E+43 .453E+31 .310E+01 .155E+31 .424E+00 .121E+02 .586E+01 .909E+00
   50 .113E+43 .356E+31 .256E+01 .122E+31 .354E+00 .114E+02 .555E+01 .915E+00
   60 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .106E+02 .520E+01 .930E+00
   70 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .998E+01 .498E+01 .945E+00
   80 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .942E+01 .476E+01 .957E+00
   90 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .888E+01 .455E+01 .966E+00
  100 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .850E+01 .482E+01 .982E+00
  120 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .768E+01 .448E+01 .994E+00
  140 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .698E+01 .417E+01 .101E+01
  160 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .638E+01 .384E+01 .101E+01
  180 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .603E+01 .364E+01 .102E+01
  200 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .571E+01 .347E+01 .103E+01
  250 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .501E+01 .307E+01 .104E+01
  300 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .442E+01 .270E+01 .105E+01
  350 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .409E+01 .252E+01 .106E+01
  400 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .374E+01 .231E+01 .107E+01
  450 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .339E+01 .210E+01 .108E+01
  500 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .320E+01 .199E+01 .109E+01
  600 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .287E+01 .181E+01 .111E+01
  700 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .257E+01 .164E+01 .112E+01
  800 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .226E+01 .146E+01 .113E+01
  900 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .209E+01 .137E+01 .114E+01
 1000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .204E+01 .136E+01 .115E+01
 1200 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .194E+01 .140E+01 .116E+01
 1400 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .152E+01 .116E+01 .116E+01
 1600 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .126E+01 .970E+00 .117E+01
 1800 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .113E+01 .872E+00 .117E+01
 2000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .105E+01 .821E+00 .118E+01
 2500 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .859E+00 .686E+00 .118E+01
 3000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .741E+00 .601E+00 .118E+01
 3500 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .662E+00 .550E+00 .118E+01
 4000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .591E+00 .501E+00 .118E+01
 4500 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .527E+00 .453E+00 .118E+01
 5000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .483E+00 .421E+00 .117E+01
 6000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .412E+00 .368E+00 .117E+01
 7000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .360E+00 .328E+00 .116E+01
 8000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .322E+00 .297E+00 .116E+01
 9000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .288E+00 .268E+00 .115E+01
10000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .267E+00 .250E+00 .115E+01
11000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .249E+00 .234E+00 .116E+01
12000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .235E+00 .220E+00 .116E+01
13000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .221E+00 .207E+00 .117E+01
14000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .209E+00 .196E+00 .117E+01
15000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .199E+00 .186E+00 .116E+01
16000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .190E+00 .178E+00 .116E+01
17000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .180E+00 .169E+00 .116E+01
18000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .170E+00 .161E+00 .116E+01
19000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .161E+00 .153E+00 .116E+01
20000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .152E+00 .146E+00 .115E+01


Ejecta of massive stars: SNII (model B of Woosley & Weaver [1995]) + stellar winds                  
Initial mass function: IMF_Salpeter.dat from   0.1000 to 100.0000 solar masses                      
Evolutionary tracks: tracksZ0.0001.dat                                                              
Metallicity (mass fraction): .0001                                                                  
Evolutionary tracks: tracksZ0.0004.dat                                                              
Metallicity (mass fraction): .0004                                                                  
Evolutionary tracks: tracksZ0.004.dat                                                               
Metallicity (mass fraction): .0040                                                                  
Evolutionary tracks: tracksZ0.008.dat                                                               
Metallicity (mass fraction): .0080                                                                  
Evolutionary tracks: tracksZ0.02.dat                                                                
Metallicity (mass fraction): .0200                                                                  
Evolutionary tracks: tracksZ0.05.dat                                                                
Metallicity (mass fraction): .0500                                                                  
Evolutionary tracks: tracksZ0.1.dat                                                                 
Metallicity (mass fraction): .1000                                                                  
Fraction of close binary systems: 0.50000E-01                                                       
Initial metallicity: 0.14200E-01                                                                    
No infall                                                                                           
Type of star formation:   0                                                                         
Consistent evolution of the stellar metallicity                                                     
Mass fraction of substellar objects: 0.00000E+00                                                    
No galactic winds                                                                                   
Nebular emission                                                                                    
No extinction                                                                                       
 time nLymcont   L(Ha)    W(Ha)    L(Hb)    W(Hb)  LB/LBsol LV/LVsol   D4000
    0 .243E+47 .337E+35 .282E+04 .116E+35 .484E+03 .151E+03 .100E+03 .686E+00
    1 .235E+47 .325E+35 .255E+04 .112E+35 .417E+03 .158E+03 .102E+03 .697E+00
    2 .193E+47 .267E+35 .193E+04 .916E+34 .285E+03 .162E+03 .972E+02 .721E+00
    3 .143E+47 .198E+35 .121E+04 .678E+34 .162E+03 .177E+03 .963E+02 .750E+00
    4 .681E+46 .944E+34 .634E+03 .324E+34 .822E+02 .145E+03 .740E+02 .790E+00
    5 .358E+46 .496E+34 .450E+03 .170E+34 .579E+02 .103E+03 .513E+02 .803E+00
    6 .106E+46 .148E+34 .918E+02 .506E+33 .140E+02 .110E+03 .587E+02 .872E+00
    7 .550E+45 .763E+33 .473E+02 .262E+33 .794E+01 .959E+02 .548E+02 .923E+00
    8 .293E+45 .406E+33 .275E+02 .139E+33 .548E+01 .687E+02 .438E+02 .945E+00
    9 .170E+45 .235E+33 .227E+02 .806E+32 .473E+01 .492E+02 .289E+02 .874E+00
   10 .105E+45 .145E+33 .169E+02 .499E+32 .341E+01 .422E+02 .246E+02 .882E+00
   12 .460E+44 .638E+32 .897E+01 .219E+32 .166E+01 .389E+02 .218E+02 .905E+00
   14 .234E+44 .324E+32 .560E+01 .111E+32 .101E+01 .328E+02 .180E+02 .901E+00
   16 .131E+44 .182E+32 .355E+01 .625E+31 .664E+00 .281E+02 .156E+02 .902E+00
   18 .802E+43 .111E+32 .239E+01 .381E+31 .457E+00 .250E+02 .139E+02 .903E+00
   20 .518E+43 .718E+31 .170E+01 .246E+31 .321E+00 .231E+02 .128E+02 .910E+00
   25 .211E+43 .293E+31 .799E+00 .100E+31 .153E+00 .198E+02 .112E+02 .927E+00
   30 .991E+42 .137E+31 .422E+00 .471E+30 .814E-01 .175E+02 .995E+01 .942E+00
   35 .527E+42 .989E+30 .347E+00 .339E+30 .658E-01 .157E+02 .890E+01 .953E+00
   40 .294E+42 .731E+30 .278E+00 .251E+30 .528E-01 .144E+02 .826E+01 .963E+00
   45 .173E+42 .544E+30 .222E+00 .186E+30 .422E-01 .135E+02 .774E+01 .971E+00
   50 .104E+42 .408E+30 .178E+00 .140E+30 .339E-01 .126E+02 .727E+01 .979E+00
   60 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .111E+02 .643E+01 .991E+00
   70 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .101E+02 .587E+01 .100E+01
   80 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .923E+01 .541E+01 .101E+01
   90 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .857E+01 .506E+01 .102E+01
  100 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .797E+01 .475E+01 .103E+01
  120 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .700E+01 .422E+01 .105E+01
  140 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .626E+01 .383E+01 .106E+01
  160 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .563E+01 .348E+01 .107E+01
  180 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .508E+01 .317E+01 .109E+01
  200 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .471E+01 .295E+01 .110E+01
  250 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .405E+01 .260E+01 .114E+01
  300 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .356E+01 .234E+01 .119E+01
  350 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .313E+01 .211E+01 .123E+01
  400 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .279E+01 .192E+01 .126E+01
  450 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .251E+01 .177E+01 .130E+01
  500 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .227E+01 .164E+01 .133E+01
  600 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .188E+01 .143E+01 .137E+01
  700 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .160E+01 .128E+01 .140E+01
  800 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .138E+01 .115E+01 .141E+01
  900 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .121E+01 .106E+01 .142E+01
 1000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .108E+01 .985E+00 .144E+01
 1200 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .896E+00 .877E+00 .145E+01
 1400 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .735E+00 .747E+00 .147E+01
 1600 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .601E+00 .629E+00 .147E+01
 1800 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .521E+00 .561E+00 .149E+01
 2000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .445E+00 .482E+00 .150E+01
 2500 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .349E+00 .396E+00 .157E+01
 3000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .292E+00 .346E+00 .163E+01
 3500 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .251E+00 .301E+00 .165E+01
 4000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .212E+00 .253E+00 .165E+01
 4500 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .191E+00 .230E+00 .167E+01
 5000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .174E+00 .212E+00 .169E+01
 6000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .150E+00 .187E+00 .174E+01
 7000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .129E+00 .162E+00 .177E+01
 8000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .115E+00 .145E+00 .180E+01
 9000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .104E+00 .134E+00 .184E+01
10000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .946E-01 .124E+00 .188E+01
11000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .870E-01 .115E+00 .191E+01
12000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .803E-01 .107E+00 .193E+01
13000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .747E-01 .100E+00 .196E+01
14000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .700E-01 .947E-01 .198E+01
15000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .658E-01 .896E-01 .200E+01
16000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .626E-01 .856E-01 .201E+01
17000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .597E-01 .819E-01 .202E+01
18000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .572E-01 .785E-01 .203E+01
19000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .550E-01 .757E-01 .204E+01
20000 .000E+00 .000E+00 .000E+00 .000E+00 .000E+00 .535E-01 .731E-01 .199E+01
#endif

#if 0
/*
 * This function loads the number of Lyman continum photon [s^-1].
 */
static void LoadnLyData(char fname[], const int step){

#define NLineHeader_nLy (26)
    FILE *fp;

    FileOpen(fp,fname,"r");
    for(int i=0;i<NLineHeader_nLy;i++){
        char log[MaxCharactersInLine];
        fgets(log,MaxCharactersInLine,fp);
    }
    int lines = 0;
    while(fscanf(fp,"%*g %*g %*g %*g %*g %*g %*g %*g %*g") != EOF){
        lines ++;
    }
    fclose(fp);

    if(step == 0)
        AgeLy = malloc(sizeof(double)*lines);
    nLy[step] = malloc(sizeof(double)*lines);


    FileOpen(fp,fname,"r");
    for(int i=0;i<NLineHeader_nLy;i++){
        char log[MaxCharactersInLine];
        fgets(log,MaxCharactersInLine,fp);
    }
    lines = 0;
    while(fscanf(fp,"%lg %lg %*g %*g %*g %*g %*g %*g %*g",AgeLy+lines,nLy[step]+lines) != EOF){
        lines ++;
    }
    fclose(fp);
    NumberofLy = lines;


    for(int i=0;i<lines;i++){
        if(nLy[step][i] == 0.e0){
            nLyEnd[step] = i;
            break;
        }
    }

    return ;
}
#endif

/*
 * This function returns the coefficient of the recombination as a function of
 * the gas temperature. The reference is the Osterbrock's book (Table 2).
 * http://astro.berkeley.edu/~ay216/08/NOTES/Lecture03-08.pdf
 * http://www.strw.leidenuniv.nl/~michiel/ismclass_files/ism08_lecture2.pdf
 */
#if 0
static double ReturnAlpha_T(const int index){
    //return (2.6e-13*pow((1e4/(Pall.ConvertUtoT*Pbody[index]->UP)),0.85));
    return (2.6e-13*pow((1e4/(Pall.ConvertUtoT*Pbody[index]->UP)),0.85));
}
#endif

double ReturnRecombinationConstantAlpha(void){
    //If we assume that the gas temperature is 1e4 K, the value is 2.6e-13
    //[s^{-1}].
    return (2.6e-13);
}

static double ReturnAlpha(void){
    //If we assume that the gas temperature is 1e4 K, the value is 2.6e-13
    //[s^{-1}].
    return (2.6e-13);
}

/*
 * This function returns the number of Lyman continum photons [s^{-1}].
 * The photon number was obtained by PEGASE.2.
 */
double ReturnNumberofLyPhoton(const int index, const int step){
    double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[index]->FormationTime)/MEGAYEAR_CGS;
    int Age_int = (int)Age;

    if(Age > AgeLy[nLyEnd[step]])
        return 0.0;

    int t_lower = 0;
    int t_upper = 0;
    for(int k=0;k<NumberofLy-1;k++){
        if(Age_int < AgeLy[k+1]){
            t_lower = k;
            t_upper = k+1;
            break;
        }
    }

    double grad = (nLy[step][t_upper]-nLy[step][t_lower])/(AgeLy[t_upper]-AgeLy[t_lower]);
    double CurrentnLy = grad*(Age - AgeLy[t_lower]) + nLy[step][t_lower];

    return CurrentnLy;
}

/*
 * This function calculates the Stromgren radius for the particle ``index''.
 * The Stromgre radius is Rs = {(3 Q)/(4 M_PI) * 1/(ne*np*\alpha) }^{1/3},
 * where Q is the number of the lyman continum photons, ne and np are the number
 * of the electron and poroton and \alpha is the recombination constant,
 * respectively. Here we assume that the gas around the source stars are fully
 * ionized for the estimation of the Stromgren radius.
 */ 
// static double ReturnStromgrenRadius(const int index){
//   
//  double Q = ReturnNumberofLyPhoton(index,0);
//  double X = 1-HeliumAbandance-Pstar[index]->Z;
//  double ne = X*Pall.ConvertNumberDensityToCGS*Pstar[index]->Density;
//  double np = ne;
//  const double fact = 3.0/(4.0*M_PI);
//
//  return cbrt(fact*Q/(ne*np*ReturnAlpha()))/Pall.UnitLength;
// }


#define KernelFactInc   (1.14) // 1.5 ^ (1.3)
#define KernelFactDec   (0.79) // 0.75 ^ (1.3)
#define KernelFactInc_First   (3) // 1.5 ^ (1.3)

struct StructHIIExport{
    double    Radius;  // Radius.
    double    Pos[3];  // Position.
    int       Leaf;
#ifdef __EXPORT_TIMESTEP__
    int       k_star;
#endif // __EXPORT_TIMESTEP__
};

struct StructHIIImport{
    double    Rho;     // Density.
#ifdef __PHOTON_COUNT_BASE__ //{
    double    PhotonCount;    // Total absorbed photon number per sec.
    double    PhotonCountDistMin;
#else // __PHOTON_COUNT_BASE__ //} //{
    double    Mass;    // Mass.
    double    MassDistMin;
#endif // __PHOTON_COUNT_BASE__ //}
    double    DistMin;
    int       Nlist;   // Nlist.
    int       Leaf;
};

struct StructActiveHIIParticle{
    int Index;
    int Nlist;
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCount;
#else //__PHOTON_COUNT_BASE__ //} //{
    double Mass;
#endif // __PHOTON_COUNT_BASE__ //}
    double LyAphoton;  // The number of Lyman continum photons [s^{-1}].
    double Radius;
    double Pos[3];

    double DistMin;
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCountDistMin;
#else //__PHOTON_COUNT_BASE__ //} //{
    double MassDistMin;
#endif // __PHOTON_COUNT_BASE__ //}
    double Rvalue;
    double Lvalue;
    bool LocalUpdateFlag;
    bool HIIRegion;
};

struct StructHIILocalInfo {
    int Nlist;
#ifdef __PHOTON_COUNT_BASE__ //{
    double PhotonCount;
    double PhotonCountDistMin;
#else //__PHOTON_COUNT_BASE__ //} //{
    double Mass;
    double MassDistMin;
#endif // __PHOTON_COUNT_BASE__ //}
    double DistMin;
};


struct StructBinarySearch {
    double Rvalue;
    double Lvalue;
};

static inline bool __attribute__((always_inline)) OverlapDomainHIIRadius(double Pos[restrict], const double h, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}


//static inline bool OverlapActiveDomainHIIradius(double Pos[restrict], const double h, const int NodeID) __attribute__((always_inline));
static inline bool __attribute__((always_inline)) OverlapActiveDomainHIIradius(double Pos[restrict], const double h, const int NodeID){ 

    double Dist2 = 0.e0;
    for(int k=0;k<3;k++){
        if(Pos[k] < EdgesForActiveHydro[NodeID].PosMin[k]) 
            Dist2 += SQ(EdgesForActiveHydro[NodeID].PosMin[k]-Pos[k]);
        if(Pos[k] > EdgesForActiveHydro[NodeID].PosMax[k])
            Dist2 += SQ(EdgesForActiveHydro[NodeID].PosMax[k]-Pos[k]);
    }
    return (Dist2 < SQ(h));
}

static struct StructHIILocalDomainEdge{
    double Max[3],Min[3];
} HIILocalDomainEdge;


static void CheckLocalDomainEdge(void){

    if(Pall.Nhydro == 0)    return;

    for(int k=0;k<3;k++){
        HIILocalDomainEdge.Max[k] = Phydro[0]->PosP[k];
        HIILocalDomainEdge.Min[k] = Phydro[0]->PosP[k];
    }
    for(int i=1;i<Pall.Nhydro;i++){
        for(int k=0;k<3;k++){
            HIILocalDomainEdge.Max[k] = fmax(Phydro[i]->PosP[k],HIILocalDomainEdge.Max[k]);
            HIILocalDomainEdge.Min[k] = fmin(Phydro[i]->PosP[k],HIILocalDomainEdge.Min[k]);
        }
    }
    return ;
}


static inline int __attribute__((always_inline)) CheckHIIExportFlags(const int Index, const int NProcs, bool HIIExportFlags[][NProcs], const int NActiveYoungStars, double Pos[restrict][3], double Raidus[restrict]){

    int ExportNodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;
    for(int i=0;i<NActiveYoungStars;i++){
        if(OverlapDomainHIIRadius(Pos[i],2.0*Raidus[i],ExportNodeID)){
            HIIExportFlags[i][Index] = ON;
            NExport ++;
        }
    }

	return NExport;
}


static inline int __attribute__((always_inline)) CheckHIIExportFlagsModified(const int NodeIndex, const int NProcs, bool HIIExportFlags[restrict], const int NActiveYoungStars, struct StructActiveHIIParticle ActiveHIIParticle[restrict]){

    if(NActiveYoungStars == 0) 
        return 0;

    int ExportNodeID = CommunicationTable[NodeIndex].SendRank;
    // Node By Node Comparison
    double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};
    if(!OverlapDomainHIIRadius(BoxCenter,LocalKernelMax,ExportNodeID)){
        return 0;
    }

    int NExport = 0;
    for(int i=0;i<NActiveYoungStars;i++){
        int Offset = i*NProcs;
        if(HIIExportFlags[Offset+NProcs-1]){
            if(OverlapDomainHIIRadius(ActiveHIIParticle[i].Pos,2.0*ActiveHIIParticle[i].Radius,ExportNodeID)){
                HIIExportFlags[Offset+NodeIndex] = true;
                NExport ++;
            }
        }
    }

	return NExport;
}


//static inline bool CheckLocalParticle(double Pos[], double Radius) __attribute__((always_inline));
static inline bool __attribute__((always_inline)) CheckLocalParticle(double Pos[], double Radius){
    int flag = true;
    for(int k=0;k<3;k++){
        if(Pos[k]+2.e0*Radius > HIILocalDomainEdge.Max[k]) flag = false;
        if(Pos[k]-2.e0*Radius < HIILocalDomainEdge.Min[k]) flag = false;
    }
    return flag;
}


/*
 * This function returns the number of neighboring particles.
 */
static int GetNumberofNeighbors(double Pos[restrict], const double Radius, int Neighbors[restrict]){

    int nlist = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,2.e0*Radius,&nlist,Neighbors);
    }while(CurrentNodeID != RootNodeID);

    return nlist;
}
#ifdef __PHOTON_COUNT_BASE__
/*
 * This function returns the absorbed photon number summing over neighboring particles.
 */

#if 0
static double ReturnNeighborsAbsorbedPhoton(double Pos[restrict], const double Radius, int Neighbors[restrict], int *nlist){

    double fact = (4.0*M_PI/3.0)*SQ(Pall.UnitMass)/CUBE(Pall.UnitLength);
    const double InvPMass = 1.0/PROTON_MASS_CGS;

    *nlist = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    double PhotonCount = 0.e0;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,2.e0*Radius,nlist,Neighbors);
        for(int i=0;i<*nlist;i++){
#ifdef __SKIP_HOTGAS__
            if(Phydro[Neighbors[i]]->UPred*Pall.ConvertUtoT > 1.e+4) continue;
#endif
            //Mass += (1.0-HeliumAbandance-Phydro[Neighbors[i]]->Z)*PhydroBody(Neighbors[i])->Mass;
            double X = (1.0-HeliumAbandance-Phydro[Neighbors[i]]->Z);
            //double Rho_H = X*Phydro[Neighbors[i]]->RhoPred/PROTON_MASS_CGS;
            double n_H = X*PhydroBody(Neighbors[i])->Mass*
                (3.0/(4.0*M_PI*CUBE(2.0*Phydro[Neighbors[i]]->KernelPred)))*InvPMass;
            double n_e = n_H;

            PhotonCount += fact*n_H*n_e*ReturnAlpha()*CUBE(2.0*Phydro[Neighbors[i]]->KernelPred);

        }
    }while(CurrentNodeID != RootNodeID);

    return PhotonCount;
}
#else
static double ReturnNeighborsAbsorbedPhoton(double Pos[restrict], const double Radius, int Neighbors[restrict], int *nlist){


    *nlist = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    double PhotonCount = 0.e0;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,2.e0*Radius,nlist,Neighbors);
        for(int i=0;i<*nlist;i++){
            int index = Neighbors[i];
#ifdef __SKIP_HOTGAS__
            if(Phydro[index]->UPred*Pall.ConvertUtoT > 1.e+4) continue;
#endif
            const double InvPMass = 1.0/PROTON_MASS_CGS;
            double X = (1.0-HeliumAbandance-Phydro[index]->Z);
            double n_H = (X*Phydro[index]->RhoPred*Pall.UnitMass/PROTON_MASS_CGS)/CUBE(Pall.UnitLength);
            /*
            double n_H = (X*PhydroBody(Neighbors[i])->Mass*Pall.UnitMass/PROTON_MASS_CGS)*3.0/
                (4.0*M_PI*CUBE(2.0*Phydro[Neighbors[i]]->KernelPred*Pall.UnitLength));
            */
            double n_e = n_H;
            double r3 = 3*PhydroBody(index)->Mass/(4*M_PI*Phydro[index]->RhoPred);
            PhotonCount += (4.0*M_PI/3.0)*n_H*n_e*ReturnAlpha()*
                r3*CUBE(Pall.UnitLength);
                //CUBE(2.0*Phydro[index]->KernelPred*Pall.UnitLength);
        }
    }while(CurrentNodeID != RootNodeID);

    return PhotonCount;
}
#endif

#else //__PHOTON_COUNT_BASE__
/*
 * This function returns the mass summing over neighboring particles.
 */
static double ReturnNeighborsMass(double Pos[restrict], const double Radius, int Neighbors[restrict], int *nlist){

    *nlist = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    double Mass = 0.e0;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,Pos,2.e0*Radius,nlist,Neighbors);
        for(int i=0;i<*nlist;i++){
#ifdef __SKIP_HOTGAS__
            if(Phydro[Neighbors[i]]->UPred*Pall.ConvertUtoT > 1.e+4) continue;
#endif
            Mass += (1.0-HeliumAbandance-Phydro[Neighbors[i]]->Z)*PhydroBody(Neighbors[i])->Mass;
        }
    }while(CurrentNodeID != RootNodeID);

    return Mass;
}
#endif // __PHOTON_COUNT_BASE__


/*
 * This function counts (1) the total number of photons which can absorbed in the neighbor particles, 
 * or (2) the total mass of the neighbor particles.
 */
static struct StructHIILocalInfo ReturnNeighborsInfo(struct StructActiveHIIParticle *ActiveHIIParticle_i, int Neighbors[restrict]){

    struct StructHIILocalInfo TempHIILocalInfo = {
        .Nlist = 0,
        .DistMin = 0.0,
        //.DistMin = 10000,
#ifdef __PHOTON_COUNT_BASE__ //{
        .PhotonCount = 0.e0,
        .PhotonCountDistMin = 0.0,
#else //__PHOTON_COUNT_BASE__ //} //{
        .Mass = 0.e0,
        .MassDistMin = 0.0,
#endif // __PHOTON_COUNT_BASE__ //}
    };

    int Iteration = 0;
    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
    do {
        CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,ActiveHIIParticle_i->Pos,
                2.e0*ActiveHIIParticle_i->Radius,&(TempHIILocalInfo.Nlist),Neighbors);


        if((TempHIILocalInfo.Nlist>0)&&(Iteration==0)){
            TempHIILocalInfo.DistMin = DISTANCE2(PhydroBody(Neighbors[0])->PosP,ActiveHIIParticle_i->Pos);
        }

        for(int i=0;i<TempHIILocalInfo.Nlist;i++){
            int index = Neighbors[i];
#ifdef __SKIP_HOTGAS__ //{
            if(Phydro[index]->UPred*Pall.ConvertUtoT > 1.e+4) continue;
#endif //__SKIP_HOTGAS__ //}

#ifdef __PHOTON_COUNT_BASE__ //{
            const double InvPMass = 1.0/PROTON_MASS_CGS;
            double X = (1.0-HeliumAbandance-Phydro[index]->Z);
            double n_H = (X*Phydro[index]->RhoPred*Pall.UnitMass/PROTON_MASS_CGS)/CUBE(Pall.UnitLength);
            //double _rho_ = 10*Phydro[index]->RhoPred/(Pall.ConvertNumberDensityToCGS*Phydro[index]->RhoPred);
            //double n_H = (X*_rho_*Pall.UnitMass/PROTON_MASS_CGS)/CUBE(Pall.UnitLength);
            double n_e = n_H;
            double r3 = 3*PhydroBody(index)->Mass/(4*M_PI*Phydro[index]->RhoPred);
            //double r3 = 3*PhydroBody(index)->Mass/(4*M_PI*_rho_);
            double CurrentPhotonCount = (4.0*M_PI/3.0)*n_H*n_e*ReturnAlpha()*r3*CUBE(Pall.UnitLength);

#if 0
            double Mass_cgs = PhydroBody(index)->Mass*Pall.UnitMass;
            double _n_e = (X*Phydro[index]->RhoPred*Pall.UnitMass/PROTON_MASS_CGS)/CUBE(Pall.UnitLength);
            double N_H = X*Mass_cgs*InvPMass*_n_e*ReturnAlpha();
            if(i == 0){
                //if(MPIGetMyID() == MPI_ROOT_RANK){
                    fprintf(stderr,"%g <> %g | %g %g %g %g %g\n",N_H,CurrentPhotonCount,
                            PhydroBody(index)->Mass,
                            PhydroBody(index)->Mass*Pall.UnitMass,
                            PhydroBody(index)->Mass*Pall.UnitMass/MSUN_CGS,
                            PROTON_MASS_CGS,X);
                //}
            }
#endif

            TempHIILocalInfo.PhotonCount += CurrentPhotonCount;
#else //__PHOTON_COUNT_BASE__ //} //{
            double CurrentMass = (1.0-HeliumAbandance-Phydro[index]->Z)*PhydroBody(index)->Mass;
            TempHIILocalInfo.Mass += CurrentMass;
#endif // __PHOTON_COUNT_BASE__ //}

            // Evaluate minimum distance.
            double TmpDistMin = DISTANCE2(PhydroBody(index)->PosP,ActiveHIIParticle_i->Pos);
            if(TempHIILocalInfo.DistMin >= TmpDistMin){
                TempHIILocalInfo.DistMin = TmpDistMin;
#ifdef __PHOTON_COUNT_BASE__ //{
                TempHIILocalInfo.PhotonCountDistMin = CurrentPhotonCount;
                //TempHIILocalInfo.PhotonCountDistMin = PhydroBody(index)->Mass*Pall.UnitMass/MSUN_CGS;
                //TempHIILocalInfo.PhotonCountDistMin = Pall.ConvertNumberDensityToCGS*Phydro[index]->RhoPred;
#else //__PHOTON_COUNT_BASE__ //} //{
                TempHIILocalInfo.MassDistMin = CurrentMass;
#endif // __PHOTON_COUNT_BASE__ //}
            }
        }
        Iteration ++;
    }while(CurrentNodeID != RootNodeID);

    return TempHIILocalInfo;
}


static inline double KernelHII(const double r, const double InvKerneli) __attribute__((always_inline));
static inline double KernelHII(const double r, const double InvKerneli){ 

	double u = r*InvKerneli;
#if (DIMENSION == 1)
    const static double coef1d = 2.0/3.0;
	double coef = coef1d*InvKerneli;
#elif (DIMENSION == 2)   
    const static double coef2d = 10.0/(7.0*M_PI);
    double coef = coef2d*SQ(InvKerneli);
#elif (DIMENSION == 3)
	double coef = M_1_PI*CUBE(InvKerneli);
#endif

	if(u<1.e0){
		return (coef*(1.e0 - 1.5*SQ(u) + 0.75*CUBE(u)));
	} else if (u<2.e0){
		return (coef*(0.25*CUBE(2.e0-u)));
	} else {
	    return 0.e0;
    }
}

static int NContactedDomains;
static int *ContactedDomainID;

static inline void AllocateContanctedDomainID(void){
    ContactedDomainID = malloc(sizeof(int)*MPIGetNumProcs());
    return ;
}

static inline bool CheckLocalExternalDomainsContacted(const int MyID, const int ExtID){

    for(int k=0;k<3;k++){
        if((EdgesForHydro[MyID].PosMax[k] < EdgesForHydro[ExtID].PosMin[k])||
           (EdgesForHydro[MyID].PosMin[k] > EdgesForHydro[ExtID].PosMax[k]))  return false;
    }
    return true;

}

/*
 * This function checkes that how many contacted domains around the local
 * domain. Comparing the local domain edge to the external domain, 
 */
static inline void CheckContactedDomain(void){
    int NProcs = MPIGetNumProcs();
    int MyID = MPIGetMyID();
    NContactedDomains = 0;
    for(int i=0;i<NProcs-1;i++){
        int NodeID = CommunicationTable[i].SendRank;
        assert(MPIGetMyID() != NodeID);
        if(CheckLocalExternalDomainsContacted(MyID,NodeID)){
            ContactedDomainID[NContactedDomains] = i;
            NContactedDomains ++;
        }
    }
    return ;
}

static inline bool CheckInLocalDomain(double Pos[], double Kernel, const int MyID) __attribute__((always_inline));
static inline bool CheckInLocalDomain(double Pos[], double Kernel, const int MyID){
    for(int k=0;k<3;k++){
        if(Pos[k]+2.e0*Kernel > EdgesForHydro[MyID].PosMax[k]) return false;
        if(Pos[k]-2.e0*Kernel < EdgesForHydro[MyID].PosMin[k]) return false;
    }
    return true;
}


#ifdef __PHOTON_COUNT_BASE__
static inline bool __attribute__((always_inline)) CheckLocalMassAndUpdateHIIradius_i(const int NProcs, bool HIIExportFlags_i[NProcs], int *NBList_i, struct StructBinarySearch *BinarySearch_i, double *Radius_i, double *PhotonCount_i, double LyAphoton_i){ 
#else // __PHOTON_COUNT_BASE__
static inline bool __attribute__((always_inline)) CheckLocalMassAndUpdateHIIradius_i(const int NProcs, bool HIIExportFlags_i[NProcs], int *NBList_i, struct StructBinarySearch *BinarySearch_i, double *Radius_i, double *Mass_i, double LyAphoton_i){ 
#endif // __PHOTON_COUNT_BASE__

#define ConvergenceFactor  (1.e-2)
    double fact = (3.0/(4.0*M_PI))*(ReturnAlpha()/(SQ(PROTON_MASS_CGS)));

#ifdef __PHOTON_COUNT_BASE__
    double Qest = *PhotonCount_i;
#else // __PHOTON_COUNT_BASE__
    double Qest = fact*(SQ((*Mass_i)*Pall.UnitMass)/(CUBE((2.0*(*Radius_i))*Pall.UnitLength)));
#endif // __PHOTON_COUNT_BASE__
    double Qmin = (1.0-ConvergenceFactor)*LyAphoton_i;
    double Qmax = (1.0+ConvergenceFactor)*LyAphoton_i;

    /*
    fprintf(stderr,"Density = %g, Ndensity = %g\n",
            *Mass_i*Pall.UnitMass/(4*M_PI/3.0)/CUBE(2.0*(*Radius_i)*Pall.UnitLength),
            *Mass_i*Pall.UnitMass/(4*M_PI/3.0)/CUBE(2.0*(*Radius_i)*Pall.UnitLength)/PROTON_MASS_CGS);
    fprintf(stderr,"R = %g, Mass = %g, np,ne = %g\n",
            2*(*Radius_i)*Pall.UnitLength/PC_CGS,*Mass_i*Pall.UnitMass/MSUN_CGS,
            *Mass_i*Pall.UnitMass/PROTON_MASS_CGS/(4*M_PI/3.0*CUBE(2.0*(*Radius_i)*Pall.UnitLength)));
    fprintf(stderr,"Qest = %g, Qmin,max = %g %g, %g %d\n",Qest,Qmin,Qmax,PhydroBody(0)->Mass,*NBList_i);
    */

    if(((Qmin)<=Qest)&&(Qest<=(Qmax))){
        HIIExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
        *PhotonCount_i = 1.0;
#else // __PHOTON_COUNT_BASE__
        *Mass_i = 1.0;
#endif // __PHOTON_COUNT_BASE__
        return true;
    }else if((Qest<=(Qmax))&&(BinarySearch_i->Rvalue>0.e0)&&(BinarySearch_i->Lvalue>0.e0)){
        if(BinarySearch_i->Rvalue-BinarySearch_i->Lvalue < 1.e-3*BinarySearch_i->Lvalue){
            HIIExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
            *PhotonCount_i = 2.0;
#else // __PHOTON_COUNT_BASE__
            *Mass_i = 2.0;
#endif // __PHOTON_COUNT_BASE__
            return true;
        }
    }
    if(HIIExportFlags_i[NProcs-1]){
        if(Qest<Qmin){
            BinarySearch_i->Lvalue = fmax(BinarySearch_i->Lvalue,*Radius_i);
        } else if(Qest>Qmax){
            if((*NBList_i) == 1){
                (*Radius_i) *= cbrt(Qmax/Qest);
                // (*Radius_i) = 0.0;
                HIIExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
                *PhotonCount_i = 3.0;
#else // __PHOTON_COUNT_BASE__
                *Mass_i = 3.0;
#endif // __PHOTON_COUNT_BASE__
                return true;
            }else{
                if(BinarySearch_i->Rvalue > 0.e0){
                    BinarySearch_i->Rvalue = fmin(BinarySearch_i->Rvalue,*Radius_i);
                }else{
                    BinarySearch_i->Rvalue = *Radius_i;
                }
            }
        }

        if((BinarySearch_i->Lvalue>0.e0)&&(BinarySearch_i->Rvalue>0.e0)){
            *Radius_i = cbrt(0.5*(CUBE(BinarySearch_i->Lvalue)+CUBE(BinarySearch_i->Rvalue)));
        }else{
            if((BinarySearch_i->Rvalue == 0.e0)&&(BinarySearch_i->Lvalue > 0.e0)){
                (*Radius_i) *= KernelFactInc_First;
            }else if((BinarySearch_i->Rvalue > 0.e0)&&(BinarySearch_i->Lvalue == 0.e0)){
                (*Radius_i) *= KernelFactDec;
            }
        }
    }
    return false;
}


#ifdef __PHOTON_COUNT_BASE__
static inline void __attribute__((always_inline)) UpdateHIILocal(double Pos_i[restrict], double *Radius_i, double *PhotonCount_i, double LyAphoton_i, int *NBList_i, int Neighbors[restrict], const int MyID, const int NProcs, bool HIIExportFlags_i[NProcs], struct StructBinarySearch *BinarySearch_i){
#else // __PHOTON_COUNT_BASE__
static inline void __attribute__((always_inline)) UpdateHIILocal(double Pos_i[restrict], double *Radius_i, double *Mass_i, double LyAphoton_i, int *NBList_i, int Neighbors[restrict], const int MyID, const int NProcs, bool HIIExportFlags_i[NProcs], struct StructBinarySearch *BinarySearch_i){
#endif //__PHOTON_COUNT_BASE__

#ifdef __PHOTON_COUNT_BASE__
    if(CheckLocalMassAndUpdateHIIradius_i(NProcs,HIIExportFlags_i,NBList_i,BinarySearch_i,Radius_i,PhotonCount_i,LyAphoton_i) == true)
        return;
#else // __PHOTON_COUNT_BASE__
    if(CheckLocalMassAndUpdateHIIradius_i(NProcs,HIIExportFlags_i,NBList_i,BinarySearch_i,Radius_i,Mass_i,LyAphoton_i) == true)
        return;
#endif //__PHOTON_COUNT_BASE__
    do{
        if(!CheckInLocalDomain(Pos_i,*Radius_i,MyID)) return;
        for(int i=0;i<NContactedDomains;i++){
            int NodeID = ContactedDomainID[i];
            if(OverlapDomainHIIRadius(Pos_i,2.0*(*Radius_i),NodeID)) return;
        }
#ifdef __PHOTON_COUNT_BASE__
        *PhotonCount_i = ReturnNeighborsAbsorbedPhoton(Pos_i,*Radius_i,Neighbors,NBList_i);
#else // __PHOTON_COUNT_BASE__
        *Mass_i = ReturnNeighborsMass(Pos_i,*Radius_i,Neighbors,NBList_i);
#endif //__PHOTON_COUNT_BASE__

#ifdef __PHOTON_COUNT_BASE__
    } while (CheckLocalMassAndUpdateHIIradius_i(NProcs,HIIExportFlags_i,NBList_i,BinarySearch_i,Radius_i,PhotonCount_i,LyAphoton_i) == false);
#else // __PHOTON_COUNT_BASE__
    } while (CheckLocalMassAndUpdateHIIradius_i(NProcs,HIIExportFlags_i,NBList_i,BinarySearch_i,Radius_i,Mass_i,LyAphoton_i) == false);
#endif //__PHOTON_COUNT_BASE__

    return;
}


static inline bool __attribute__((always_inline)) CheckLocalMassAndUpdateHIIRadiusModified_i(const int NProcs, bool HIIExportFlags_i[NProcs], struct StructActiveHIIParticle *ActiveHIIParticle_i){ 

#define ConvergenceFactor  (1.e-2)

#ifdef __PHOTON_COUNT_BASE__
    double Qest = ActiveHIIParticle_i->PhotonCount;
#else // __PHOTON_COUNT_BASE__
    double fact = (3.0/(4.0*M_PI))*(ReturnAlpha()/(SQ(PROTON_MASS_CGS)));
    double Qest = fact*(SQ((ActiveHIIParticle_i->Mass)*Pall.UnitMass)/(CUBE((2.0*(ActiveHIIParticle_i->Radius))*Pall.UnitLength)));
#endif // __PHOTON_COUNT_BASE__
    double Qmin = (1.0-ConvergenceFactor)*ActiveHIIParticle_i->LyAphoton;
    double Qmax = (1.0+ConvergenceFactor)*ActiveHIIParticle_i->LyAphoton;


    if(((Qmin)<=Qest)&&(Qest<=(Qmax))){
        HIIExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
        ActiveHIIParticle_i->PhotonCount = 1.0;
#else // __PHOTON_COUNT_BASE__
        ActiveHIIParticle_i->Mass = 1.0;
#endif // __PHOTON_COUNT_BASE__
        return true;
    }else if((Qest<=(Qmax))&&(ActiveHIIParticle_i->Rvalue>0.e0)&&(ActiveHIIParticle_i->Lvalue>0.e0)){
        if(ActiveHIIParticle_i->Rvalue-ActiveHIIParticle_i->Lvalue < 1.e-3*ActiveHIIParticle_i->Lvalue){
            HIIExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
            ActiveHIIParticle_i->PhotonCount = 2.0;
#else // __PHOTON_COUNT_BASE__ // { //}
            ActiveHIIParticle_i->Mass = 2.0;
#endif // __PHOTON_COUNT_BASE__ //}
            return true;
        }
    } 
#ifdef __PHOTON_COUNT_BASE__ //{
        else if((ActiveHIIParticle_i->PhotonCountDistMin > ActiveHIIParticle_i->LyAphoton))
#else // __PHOTON_COUNT_BASE__ // { //}
#error USE PHOTON COUNT BASE MODE
#endif // __PHOTON_COUNT_BASE__ //}
        {
            HIIExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
            ActiveHIIParticle_i->PhotonCount = 4.0;
#else // __PHOTON_COUNT_BASE__ // { //}
            ActiveHIIParticle_i->Mass = 4.0;
#endif // __PHOTON_COUNT_BASE__ //}
            ActiveHIIParticle_i->HIIRegion = false;
            return true;
    }


    if(HIIExportFlags_i[NProcs-1]){
        if(Qest<Qmin){
            ActiveHIIParticle_i->Lvalue = fmax(ActiveHIIParticle_i->Lvalue,ActiveHIIParticle_i->Radius);
        } else if(Qest>Qmax){
            if((ActiveHIIParticle_i->Nlist) == 1){
                ActiveHIIParticle_i->Radius *= cbrt(Qmax/Qest);
                HIIExportFlags_i[NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
                ActiveHIIParticle_i->PhotonCount = 3.0;
#else // __PHOTON_COUNT_BASE__
                ActiveHIIParticle_i->Mass = 3.0;
#endif // __PHOTON_COUNT_BASE__
                return true;
            }else{
                if(ActiveHIIParticle_i->Rvalue > 0.e0){
                    ActiveHIIParticle_i->Rvalue = fmin(ActiveHIIParticle_i->Rvalue,ActiveHIIParticle_i->Radius);
                }else{
                    ActiveHIIParticle_i->Rvalue = ActiveHIIParticle_i->Radius;
                }
            }
        }

        if((ActiveHIIParticle_i->Lvalue>0.e0)&&(ActiveHIIParticle_i->Rvalue>0.e0)){
            ActiveHIIParticle_i->Radius = cbrt(0.5*(CUBE(ActiveHIIParticle_i->Lvalue)+CUBE(ActiveHIIParticle_i->Rvalue)));
        }else{
            if((ActiveHIIParticle_i->Rvalue == 0.e0)&&(ActiveHIIParticle_i->Lvalue > 0.e0)){
                ActiveHIIParticle_i->Radius *= KernelFactInc_First;
            }else if((ActiveHIIParticle_i->Rvalue > 0.e0)&&(ActiveHIIParticle_i->Lvalue == 0.e0)){
                ActiveHIIParticle_i->Radius *= KernelFactDec;
            }
        }
    }
    return false;
}

static inline void __attribute__((always_inline)) UpdateHIIRadiusLocalModified(const int Index, const int NActives, 
        struct StructActiveHIIParticle ActiveHIIParticle[restrict], int Neighbors[restrict], const int MyID, const int NProcs, 
            bool HIIExportFlags[restrict]){

    if(CheckLocalMassAndUpdateHIIRadiusModified_i(NProcs,HIIExportFlags+Index*NProcs,ActiveHIIParticle+Index) == true)
        return ;

    do{
        if(!CheckInLocalDomain(ActiveHIIParticle[Index].Pos,ActiveHIIParticle[Index].Radius,MyID)) return;
        for(int i=0;i<NContactedDomains;i++){
            int NodeID = ContactedDomainID[i];
            if(OverlapDomainHIIRadius(ActiveHIIParticle[Index].Pos,2.0*ActiveHIIParticle[Index].Radius,NodeID)) return;
        }

        struct StructHIILocalInfo TemporalData =
            ReturnNeighborsInfo(ActiveHIIParticle+Index,Neighbors);

            ActiveHIIParticle[Index].Nlist = TemporalData.Nlist;
#ifdef __PHOTON_COUNT_BASE__
            ActiveHIIParticle[Index].PhotonCount = TemporalData.PhotonCount;
#else //__PHOTON_COUNT_BASE__
            ActiveHIIParticle[Index].Mass = TemporalData.Mass;
#endif // __PHOTON_COUNT_BASE__
            if(TemporalData.Nlist > 0){
                ActiveHIIParticle[Index].DistMin = TemporalData.DistMin;
#ifdef __PHOTON_COUNT_BASE__
                ActiveHIIParticle[Index].PhotonCountDistMin = TemporalData.PhotonCountDistMin;
#else //__PHOTON_COUNT_BASE__
                ActiveHIIParticle[Index].MassDistMin = TemporalData.MassDistMin;
#endif // __PHOTON_COUNT_BASE__
            }
    } while (CheckLocalMassAndUpdateHIIRadiusModified_i(NProcs,HIIExportFlags+Index*NProcs,ActiveHIIParticle+Index) == false);

    return;
}


static inline int CheckNeighborNumberAndUpdateHIIradius(const int NActives, const int NProcs, bool HIIExportFlags[restrict][NProcs], int ActiveIndexList[restrict], int ActiveNeighborList[restrict], struct StructBinarySearch BinarySearch[restrict], double Radius[restrict])  __attribute__((always_inline));
static inline int CheckNeighborNumberAndUpdateHIIradius(const int NActives, const int NProcs, bool HIIExportFlags[restrict][NProcs], int ActiveIndexList[restrict], int ActiveNeighborList[restrict], struct StructBinarySearch BinarySearch[restrict], double Radius[]){ 

    int NBmin = Pall.Ns-Pall.Npm;
    int NBmax = Pall.Ns+Pall.Npm;

    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        if(HIIExportFlags[i][NProcs-1]){ 
            int Nlist = ActiveNeighborList[i];
            if(((NBmin)<=Nlist)&&(Nlist<=(NBmax))){
                HIIExportFlags[i][NProcs-1] = false;
            }else if((Nlist<=(NBmax))&&(BinarySearch[i].Rvalue>0.e0)&&(BinarySearch[i].Lvalue>0.e0)){
                if(BinarySearch[i].Rvalue-BinarySearch[i].Lvalue < 1.e-3*BinarySearch[i].Lvalue)
                    HIIExportFlags[i][NProcs-1] = false;
            }
            if(HIIExportFlags[i][NProcs-1]){
                if(Nlist<NBmin){
                    BinarySearch[i].Lvalue = fmax(BinarySearch[i].Lvalue,Radius[i]);
                } else if(Nlist>NBmax){
                    if(BinarySearch[i].Rvalue > 0.e0){
                        BinarySearch[i].Rvalue = fmin(BinarySearch[i].Rvalue,Radius[i]);
                    }else{
                        BinarySearch[i].Rvalue = Radius[i];
                    }
                }

                if((BinarySearch[i].Lvalue>0.e0)&&(BinarySearch[i].Rvalue>0.e0)){
                    Radius[i] = cbrt(0.5*(CUBE(BinarySearch[i].Lvalue)+CUBE(BinarySearch[i].Rvalue)));
                }else{
                    if((BinarySearch[i].Rvalue == 0.e0)&&(BinarySearch[i].Lvalue > 0.e0)){
                        Radius[i] *= KernelFactInc;
                    }else if((BinarySearch[i].Rvalue > 0.e0)&&(BinarySearch[i].Lvalue == 0.e0)){
                        Radius[i] *= KernelFactDec;
                    }
                }
                NLocalActiveLeaves ++;
            }
        }
    }
    /*
    for(int i=0;i<NActives;i++){
        if(HIIExportFlags[i][NProcs-1]){ 
            int leaf = ActiveIndexList[i];
            Pstar[leaf]->StromgrenRadius = Radius[i];
        }
    }*/

    return NLocalActiveLeaves;
}

#ifdef __PHOTON_COUNT_BASE__
static inline int __attribute__((always_inline)) CheckLocalMassAndUpdateHIIradius(const int NActives, const int NProcs, bool HIIExportFlags[restrict][NProcs], int ActiveIndexList[restrict], int ActiveNeighborList[restrict], struct StructBinarySearch BinarySearch[restrict], double Radius[restrict], double PhotonCount[restrict], double LyAphoton[restrict], bool LocalUpdateFlags[restrict]){ 
#else // __PHOTON_COUNT_BASE__
static inline int __attribute__((always_inline)) CheckLocalMassAndUpdateHIIradius(const int NActives, const int NProcs, bool HIIExportFlags[restrict][NProcs], int ActiveIndexList[restrict], int ActiveNeighborList[restrict], struct StructBinarySearch BinarySearch[restrict], double Radius[restrict], double Mass[restrict], double LyAphoton[restrict], bool LocalUpdateFlags[restrict]){ 
#endif // __PHOTON_COUNT_BASE__

#define ConvergenceFactor  (1.e-2)

    // If the mass of an ionized region is smaller than that of the nearest SPH
    // particle, a probabilistic method is used to determine the ionized region.

    double fact = (3.0/(4.0*M_PI))*(ReturnAlpha()/(SQ(PROTON_MASS_CGS)));
    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        if(HIIExportFlags[i][NProcs-1]){ 
#ifdef __LOCAL__
        if(LocalUpdateFlags[i] == false){ 
#endif
#ifdef __PHOTON_COUNT_BASE__
            double Qest = PhotonCount[i];
#else // __PHOTON_COUNT_BASE__
            double Qest = fact*(SQ(Mass[i]*Pall.UnitMass)/(CUBE(2.0*Radius[i]*Pall.UnitLength)));
#endif // __PHOTON_COUNT_BASE__

            double Qmin = (1.0-ConvergenceFactor)*LyAphoton[i];
            double Qmax = (1.0+ConvergenceFactor)*LyAphoton[i];

            if(((Qmin)<=Qest)&&(Qest<=(Qmax))){
                HIIExportFlags[i][NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
                PhotonCount[i] = 1.0;
#else // __PHOTON_COUNT_BASE__
                Mass[i] = 1.0;
#endif // __PHOTON_COUNT_BASE__
            }else if((Qest<=(Qmax))&&(BinarySearch[i].Rvalue>0.e0)&&(BinarySearch[i].Lvalue>0.e0)){
                if(BinarySearch[i].Rvalue-BinarySearch[i].Lvalue < 1.e-3*BinarySearch[i].Lvalue){
                    HIIExportFlags[i][NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
                    PhotonCount[i] = 2.0;
#else // __PHOTON_COUNT_BASE__
                    Mass[i] = 2.0;
#endif // __PHOTON_COUNT_BASE__
                }
            }
            if(HIIExportFlags[i][NProcs-1]){
                if(Qest<Qmin){
                    BinarySearch[i].Lvalue = fmax(BinarySearch[i].Lvalue,Radius[i]);
                } else if(Qest>Qmax){
                    if(ActiveNeighborList[i] == 1){
                    // if((ActiveNeighborList[i] < 32)&&(ActiveNeighborList[i]>0)){
                        //x = Rs*cbrt(Qmax/Qest);
                        Radius[i] *= cbrt(Qmax/Qest);
                        // Radius[i] = 0.0;
                        HIIExportFlags[i][NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__
                        PhotonCount[i] = 3.0;
#else // __PHOTON_COUNT_BASE__
                        Mass[i] = 3.0;
#endif // __PHOTON_COUNT_BASE__
                    }else{
                        if(BinarySearch[i].Rvalue > 0.e0){
                            BinarySearch[i].Rvalue = fmin(BinarySearch[i].Rvalue,Radius[i]);
                        }else{
                            BinarySearch[i].Rvalue = Radius[i];
                        }
                    }
                }

                if(HIIExportFlags[i][NProcs-1]){
                if((BinarySearch[i].Lvalue>0.e0)&&(BinarySearch[i].Rvalue>0.e0)){
                    Radius[i] = cbrt(0.5*(CUBE(BinarySearch[i].Lvalue)+CUBE(BinarySearch[i].Rvalue)));
                }else{
                    if((BinarySearch[i].Rvalue == 0.e0)&&(BinarySearch[i].Lvalue > 0.e0)){
                        //Radius[i] *= KernelFactInc;
                        Radius[i] *= KernelFactInc_First;
                    }else if((BinarySearch[i].Rvalue > 0.e0)&&(BinarySearch[i].Lvalue == 0.e0)){
                        Radius[i] *= KernelFactDec;
                    }
                }
                NLocalActiveLeaves ++;
                }
            }
#ifdef __LOCAL__
        } else {
            NLocalActiveLeaves ++;
        }
#endif
        }
    }

    return NLocalActiveLeaves;
}

#ifdef __PHOTON_COUNT_BASE__
static inline int __attribute__((always_inline)) CheckLocalMassAndUpdateHIIradiusModified(const int NActives, const int NProcs, bool HIIExportFlags[restrict], int ActiveIndexList[restrict], int ActiveNeighborList[restrict], struct StructBinarySearch BinarySearch[restrict], double Radius[restrict], double PhotonCount[restrict], double LyAphoton[restrict], bool LocalUpdateFlags[restrict]){ 
#else // __PHOTON_COUNT_BASE__
static inline int __attribute__((always_inline)) CheckLocalMassAndUpdateHIIradiusModified(const int NActives, const int NProcs, bool HIIExportFlags[restrict], int ActiveIndexList[restrict], int ActiveNeighborList[restrict], struct StructBinarySearch BinarySearch[restrict], double Radius[restrict], double Mass[restrict], double LyAphoton[restrict], bool LocalUpdateFlags[restrict]){ 
#endif // __PHOTON_COUNT_BASE__

#define ConvergenceFactor  (1.e-2)

    // If the mass of an ionized region is smaller than that of the nearest SPH
    // particle, a probabilistic method is used to determine the ionized region.

    double fact = (3.0/(4.0*M_PI))*(ReturnAlpha()/(SQ(PROTON_MASS_CGS)));
    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(HIIExportFlags[Offset+NProcs-1]){ 
#ifdef __LOCAL__ //{
        if(LocalUpdateFlags[i] == false){ 
#endif //__LOCAL__ //}

#ifdef __PHOTON_COUNT_BASE__ //{
            double Qest = PhotonCount[i];
#else // __PHOTON_COUNT_BASE__ //{ //}
            double Qest = fact*(SQ(Mass[i]*Pall.UnitMass)/(CUBE(2.0*Radius[i]*Pall.UnitLength)));
#endif // __PHOTON_COUNT_BASE__ //}

            double Qmin = (1.0-ConvergenceFactor)*LyAphoton[i];
            double Qmax = (1.0+ConvergenceFactor)*LyAphoton[i];

            if(((Qmin)<=Qest)&&(Qest<=(Qmax))){
                HIIExportFlags[Offset+NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
                PhotonCount[i] = 1.0;
#else // __PHOTON_COUNT_BASE__ //} //{
                Mass[i] = 1.0;
#endif // __PHOTON_COUNT_BASE__ //}
            }else if((Qest<=(Qmax))&&(BinarySearch[i].Rvalue>0.e0)&&(BinarySearch[i].Lvalue>0.e0)){
                if(BinarySearch[i].Rvalue-BinarySearch[i].Lvalue < 1.e-3*BinarySearch[i].Lvalue){
                    HIIExportFlags[Offset+NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
                    PhotonCount[i] = 2.0;
#else // __PHOTON_COUNT_BASE__ //{ //}
                    Mass[i] = 2.0;
#endif // __PHOTON_COUNT_BASE__ //}
                }
            }

            if(HIIExportFlags[Offset+NProcs-1]){
                if(Qest<Qmin){
                    BinarySearch[i].Lvalue = fmax(BinarySearch[i].Lvalue,Radius[i]);
                } else if(Qest>Qmax){
                    if(ActiveNeighborList[i] == 1){
                    // if((ActiveNeighborList[i] < 32)&&(ActiveNeighborList[i]>0)){
                        //x = Rs*cbrt(Qmax/Qest);
                        Radius[i] *= cbrt(Qmax/Qest);
                        // Radius[i] = 0.0;
                        //HIIExportFlags[i][NProcs-1] = false;
                        HIIExportFlags[Offset+NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
                        PhotonCount[i] = 3.0;
#else // __PHOTON_COUNT_BASE__ //} //{
                        Mass[i] = 3.0;
#endif // __PHOTON_COUNT_BASE__ //}
                    }else{
                        if(BinarySearch[i].Rvalue > 0.e0){
                            BinarySearch[i].Rvalue = fmin(BinarySearch[i].Rvalue,Radius[i]);
                        }else{
                            BinarySearch[i].Rvalue = Radius[i];
                        }
                    }
                }

                if(HIIExportFlags[Offset+NProcs-1]){
                    if((BinarySearch[i].Lvalue>0.e0)&&(BinarySearch[i].Rvalue>0.e0)){
                        Radius[i] = cbrt(0.5*(CUBE(BinarySearch[i].Lvalue)+CUBE(BinarySearch[i].Rvalue)));
                    }else{
                        if((BinarySearch[i].Rvalue == 0.e0)&&(BinarySearch[i].Lvalue > 0.e0)){
                            Radius[i] *= KernelFactInc_First;
                        }else if((BinarySearch[i].Rvalue > 0.e0)&&(BinarySearch[i].Lvalue == 0.e0)){
                            Radius[i] *= KernelFactDec;
                        }
                    }
                    NLocalActiveLeaves ++;
                }
            }

            // fprintf(stderr,"Check HII %d %ld %d %g | %g %g %g | %g\n",ActiveIndexList[i],PstarBody(ActiveIndexList[i])->GlobalID,
                    // ActiveNeighborList[i],Radius[i],Qmin,Qest,Qmax,PhotonCount[i]);
            // fflush(NULL);

#ifdef __LOCAL__ //{
        } else {
            NLocalActiveLeaves ++;
        }
#endif //__LOCAL__ //}
        }
    }

    return NLocalActiveLeaves;
}

static inline int __attribute__((always_inline)) CheckLocalMassAndUpdateHIIRadiusModified2(const int NActives, const int NProcs, bool HIIExportFlags[restrict], struct StructActiveHIIParticle ActiveHIIParticle[restrict]){ 

#define ConvergenceFactor  (1.e-2)

    // If the mass of an ionized region is smaller than that of the nearest SPH
    // particle, a probabilistic method is used to determine the ionized region.

    double fact = (3.0/(4.0*M_PI))*(ReturnAlpha()/(SQ(PROTON_MASS_CGS)));
    int NLocalActiveLeaves = 0;
    for(int i=0;i<NActives;i++){
        int Offset = i*NProcs;
        if(HIIExportFlags[Offset+NProcs-1]){ 
#ifdef __LOCAL__ //{
        if(ActiveHIIParticle[i].LocalUpdateFlag == false){ 
#endif //__LOCAL__ //}

#ifdef __PHOTON_COUNT_BASE__ //{
            double Qest = ActiveHIIParticle[i].PhotonCount;
#else // __PHOTON_COUNT_BASE__ //{ //}
            double Qest = fact*(SQ(ActiveHIIParticle[i].Mass*Pall.UnitMass)/(CUBE(2.0*ActiveHIIParticle[i].Radius*Pall.UnitLength)));
#endif // __PHOTON_COUNT_BASE__ //}

            double Qmin = (1.0-ConvergenceFactor)*ActiveHIIParticle[i].LyAphoton;
            double Qmax = (1.0+ConvergenceFactor)*ActiveHIIParticle[i].LyAphoton;

            if(((Qmin)<=Qest)&&(Qest<=(Qmax))){
                HIIExportFlags[Offset+NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
                ActiveHIIParticle[i].PhotonCount = 1.0;
#else // __PHOTON_COUNT_BASE__ //} //{
                ActiveHIIParticle[i].Mass = 1.0;
#endif // __PHOTON_COUNT_BASE__ //}
            }else if((Qest<=(Qmax))&&(ActiveHIIParticle[i].Rvalue>0.e0)&&(ActiveHIIParticle[i].Lvalue>0.e0)){
                if(ActiveHIIParticle[i].Rvalue-ActiveHIIParticle[i].Lvalue < 1.e-3*ActiveHIIParticle[i].Lvalue){
                    HIIExportFlags[Offset+NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
                    ActiveHIIParticle[i].PhotonCount = 2.0;
#else // __PHOTON_COUNT_BASE__ //{ //}
                    ActiveHIIParticle[i].Mass = 2.0;
#endif // __PHOTON_COUNT_BASE__ //}
                }
            }
#ifdef __PHOTON_COUNT_BASE__ //{
                else if((ActiveHIIParticle[i].PhotonCountDistMin > ActiveHIIParticle[i].LyAphoton))
#else // __PHOTON_COUNT_BASE__ // { //}
#endif // __PHOTON_COUNT_BASE__ //}
                {
                    HIIExportFlags[Offset+NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
                    ActiveHIIParticle[i].PhotonCount = 4.0;
#else // __PHOTON_COUNT_BASE__ // { //}
                    ActiveHIIParticle[i].Mass = 4.0;
#endif // __PHOTON_COUNT_BASE__ //}
                    ActiveHIIParticle[i].HIIRegion = false;
            }

            if(HIIExportFlags[Offset+NProcs-1]){
                if(Qest<Qmin){
                    ActiveHIIParticle[i].Lvalue = fmax(ActiveHIIParticle[i].Lvalue,ActiveHIIParticle[i].Radius);
                } else if(Qest>Qmax){
                    if(ActiveHIIParticle[i].Nlist == 1){
                        ActiveHIIParticle[i].Radius *= cbrt(Qmax/Qest);
                        HIIExportFlags[Offset+NProcs-1] = false;
#ifdef __PHOTON_COUNT_BASE__ //{
                        ActiveHIIParticle[i].PhotonCount = 3.0;
#else // __PHOTON_COUNT_BASE__ //} //{
                        ActiveHIIParticle[i].Mass = 3.0;
#endif // __PHOTON_COUNT_BASE__ //}
                    }else{
                        if(ActiveHIIParticle[i].Rvalue > 0.e0){
                            ActiveHIIParticle[i].Rvalue = fmin(ActiveHIIParticle[i].Rvalue,ActiveHIIParticle[i].Radius);
                        }else{
                            ActiveHIIParticle[i].Rvalue = ActiveHIIParticle[i].Radius;
                        }
                    }
                }

                if(HIIExportFlags[Offset+NProcs-1]){
                    if((ActiveHIIParticle[i].Lvalue>0.e0)&&(ActiveHIIParticle[i].Rvalue>0.e0)){
                        ActiveHIIParticle[i].Radius = cbrt(0.5*(CUBE(ActiveHIIParticle[i].Lvalue)+CUBE(ActiveHIIParticle[i].Rvalue)));
                    }else{
                        if((ActiveHIIParticle[i].Rvalue == 0.e0)&&(ActiveHIIParticle[i].Lvalue > 0.e0)){
                            ActiveHIIParticle[i].Radius *= KernelFactInc_First;
                        }else if((ActiveHIIParticle[i].Rvalue > 0.e0)&&(ActiveHIIParticle[i].Lvalue == 0.e0)){
                            ActiveHIIParticle[i].Radius *= KernelFactDec;
                        }
                    }
                    NLocalActiveLeaves ++;
                }
            }

            /*
            fprintf(stderr,"Check HII %d %ld %d %g | %g %g %g | %g\n",ActiveIndexList[i],PstarBody(ActiveIndexList[i])->GlobalID,
                    ActiveNeighborList[i],Radius[i],Qmin,Qest,Qmax,PhotonCount[i]);
            fflush(NULL);
            */

#ifdef __LOCAL__ //{
        } else {
            NLocalActiveLeaves ++;
        }
#endif //__LOCAL__ //}
        }
    }

    return NLocalActiveLeaves;
}

int CheckHIIflags(const int step){

    int NLocalActiveYoungStars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        Pstar[i]->HIIflag = false;
        if(!Pstar[i]->TypeII){
#ifdef PRESERVE_SNII_EVENTRATE //{
        if(Pstar[i]->TypeIIProb){
            double MassInMsun = Pstar[i]->InitialMass*Pall.UnitMass/MSUN_CGS;
            double prob = SNIINumber*MassInMsun;
            if(!((prob >= 1.0) || (Pstar[i]->TypeIIProb)))
                continue;
#endif // PRESERVE_SNII_EVENTRATE //}
#ifdef __UPDATE_ONLY_ACTIVE__ //{
        if((PstarBody(i)->Active == true) || (Pstar[i]->StromgrenRadius == 0.e0)){
#endif // __UPDATE_ONLY_ACTIVE__ //}//{
            double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[i]->FormationTime)/MEGAYEAR_CGS;
            if(Age < AgeLy[nLyEnd[step]]){
                Pstar[i]->HIIflag = true;
                NLocalActiveYoungStars ++;
            }
#ifdef __UPDATE_ONLY_ACTIVE__ //{
        }
#endif //__UPDATE_ONLY_ACTIVE__ //}
#ifdef PRESERVE_SNII_EVENTRATE //{
        }
#endif //}
        }
    }

    return NLocalActiveYoungStars;
}


static void WriteHIIradius(void){
    FILE *fp;
    FileOpen(fp,"./HIIradius.dat","a");

    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->HIIflag){
            fprintf(fp,"%d %ld %g %g %g %g %g\n",i, PstarBody(i)->GlobalID, 
                    Pstar[i]->StromgrenRadius*Pall.UnitLength/PC_CGS,
                    //Pstar[i]->Density*Pall.ConvertNumberDensityToCGS,
                    Pstar[i]->Density, // This is the HII region check condition.
                    PstarBody(i)->PosP[0],PstarBody(i)->PosP[1],PstarBody(i)->PosP[2]);
        }
    }
    fflush(fp);
    fclose(fp);

    return ;
}


static inline int CheckActiveHIIExportFlags(const int Index, const int NProcs, bool HIIExportFlags[restrict][NProcs], const int NActiveYoungStars, double Pos[restrict][3], double Raidus[restrict]) __attribute__((always_inline));
static inline int CheckActiveHIIExportFlags(const int Index, const int NProcs, bool HIIExportFlags[restrict][NProcs], const int NActiveYoungStars, double Pos[restrict][3], double Raidus[restrict]){

    //if(Pall.Nhydro == 0)
        //return 0;

    int ExportNodeID = CommunicationTable[Index].SendRank;
    int NExport = 0;
    for(int i=0;i<NActiveYoungStars;i++){
        if(OverlapActiveDomainHIIradius(Pos[i],2.0*Raidus[i],ExportNodeID)){
            HIIExportFlags[i][Index] = true;
            NExport ++;
        }
    }

	return NExport;
}


/*
 * When a gas particle is in the Stromgren radii of young stars, turn on the
 * HIIflag of the gas particle.
 */
static void MakeHIIregion(void){

    int NProcs = MPIGetNumProcs();
    int Neighbors[MaxNeighborSize];
    MPI_Status  mpi_status;

    int NActiveYoungStars = 0;
#ifdef __UPDATE_ONLY_ACTIVE__
    for(int i=0;i<Pall.Nstars;i++){
        Pstar[i]->HIIflag = false;
#ifdef PRESERVE_SNII_EVENTRATE
        if((!Pstar[i]->TypeII)&&(Pstar[i]->TypeIIProb)){
#else
        if(!Pstar[i]->TypeII){
#endif
            if(Pstar[i]->StromgrenRadius < 0.0) continue;

            double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[i]->FormationTime)/MEGAYEAR_CGS;
            if(Age < AgeLy[nLyEnd[0]]){
                Pstar[i]->HIIflag = true;
                NActiveYoungStars ++;
            }
        }
    }
#else
    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->HIIflag)
            NActiveYoungStars ++;
    }
#endif



    // Get the list of the young stars from external domains.
    static int HIIExportFlagsMaxAllocated = 0;
    static bool (*HIIExportFlags)[NProcs];
    static double (*Pos)[3];
    static double *Radius;
#ifdef __EXPORT_TIMESTEP__
    static int *k_star;
#endif //__EXPORT_TIMESTEP__
    if(HIIExportFlagsMaxAllocated < MAX(ForAngelsShare*NActiveYoungStars,NAdditionUnit)){
        if(HIIExportFlagsMaxAllocated > 0){
            free(HIIExportFlags);
            free(Radius);
            free(Pos);
#ifdef __EXPORT_TIMESTEP__
            free(k_star);
#endif //__EXPORT_TIMESTEP__
        }
        HIIExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*NActiveYoungStars,NAdditionUnit));
        HIIExportFlags = malloc(sizeof(bool)*HIIExportFlagsMaxAllocated*NProcs);
        Radius = malloc(sizeof(double)*HIIExportFlagsMaxAllocated);
        Pos = malloc(sizeof(double)*HIIExportFlagsMaxAllocated*3);
#ifdef __EXPORT_TIMESTEP__
        k_star = malloc(sizeof(int)*HIIExportFlagsMaxAllocated);
#endif //__EXPORT_TIMESTEP__
    }

    NActiveYoungStars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->HIIflag){
            Pos[NActiveYoungStars][0] = PstarBody(i)->PosP[0];
            Pos[NActiveYoungStars][1] = PstarBody(i)->PosP[1];
            Pos[NActiveYoungStars][2] = PstarBody(i)->PosP[2];
            Radius[NActiveYoungStars] = Pstar[i]->StromgrenRadius;
            //assert(Radius[NActiveYoungStars]<10*PC_CGS/Pall.UnitLength);
#ifdef __EXPORT_TIMESTEP__
            // k_star[NActiveYoungStars] = PstarBody(i)->k-__DT_DIFF__;
            k_star[NActiveYoungStars] = PstarBody(i)->k;
#endif //__EXPORT_TIMESTEP__
                
            for(int k=0;k<NProcs;k++)
                HIIExportFlags[NActiveYoungStars][k] = false;
            NActiveYoungStars ++;
        }
    }

    int NExportThisTime[NProcs];
    int NImportThisTime[NProcs];

    struct StructHIIExport *HIIExportSend[NProcs];
    struct StructHIIExport *HIIExportRecv = NULL;
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    for(int i=0;i<NProcs-1;i++){
        NExportThisTime[i] = CheckActiveHIIExportFlags(i,NProcs,HIIExportFlags,NActiveYoungStars,Pos,Radius);
        //fprintf(stderr,"[%d] %d\n",MPIGetMyID(),NExportThisTime[i]);
        //fflush(stderr);
        //CheckSizeofBufferExportSendIndex(_NExportThisTime[i],sizeof(struct StructHIIExport),i);
        CheckSizeofBufferExportSendIndex(NActiveYoungStars,sizeof(struct StructHIIExport),i);
        HIIExportSend[i] = BufferExportSend[i];

        int NExport = 0;
        if(NExportThisTime[i]>0){
        for(int k=0;k<NActiveYoungStars;k++){
            if(HIIExportFlags[k][i]){ 
                HIIExportSend[i][NExport].Pos[0] = Pos[k][0];
                HIIExportSend[i][NExport].Pos[1] = Pos[k][1];
                HIIExportSend[i][NExport].Pos[2] = Pos[k][2];
                HIIExportSend[i][NExport].Radius = Radius[k];
#ifdef __EXPORT_TIMESTEP__ //{
                HIIExportSend[i][NExport].k_star = k_star[k];
                //HIIExportSend[i][NExport].Leaf = k;
#endif // __EXPORT_TIMESTEP__ //}
                NExport ++;
            }
        }
        }
        assert(NExport == NExportThisTime[i]);
        NExportThisTime[i] = NExport;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    int NImport = 0;
    int NImportThisTime2[NProcs];
    int NExportThisTime2[NProcs];
    NImportThisTime2[MPIGetMyID()] = 0;
    for(int i=0;i<NProcs-1;i++){
        NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
    }
    MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
    for(int i=0;i<NProcs-1;i++){
        NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
        NImport += NImportThisTime[i];
    }
    int NImportAll = NImport;
    
    CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructHIIExport));
    HIIExportRecv = BufferExportRecv;

    NImport = 0;
    int counter_send = 0;
    int counter_recv = 0;

    int SendFlag,RecvFlag;
    for(int i=0;i<NProcs-1;i++){
        if(NExportThisTime[i]>0){
            MPI_Isend(HIIExportSend[i],
                NExportThisTime[i]*sizeof(struct StructHIIExport),
                    MPI_BYTE,CommunicationTable[i].SendRank,TAG_HII_DENSITY_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
            MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
            counter_send ++;
        }
        if(NImportThisTime[i]>0){
            MPI_Irecv(HIIExportRecv+NImport,
                NImportThisTime[i]*sizeof(struct StructHIIExport),
                    MPI_BYTE,CommunicationTable[i].RecvRank,TAG_HII_DENSITY_EXPORT+i,
                        MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
            MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
            counter_recv ++;
            NImport += NImportThisTime[i];
        }
    }
    
#ifdef TASK_TEST_STROMGRENSPHERE
    for(int i=0;i<Pall.Nhydro;i++)
        HIIFlag[i] = false;
#endif //TASK_TEST_STROMGRENSPHERE
 
    for(int i=0;i<Pall.Nstars;i++){ // for local
        if(Pstar[i]->HIIflag){
            int nlist = 0;
            int RootNodeID = 0;
            int CurrentNodeID = HydroNode[RootNodeID].Children;
            do {
                CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,
                        PstarBody(i)->PosP,Pstar[i]->StromgrenRadius,&nlist,Neighbors);

                for(int k=0;k<nlist;k++){
                    int index = Neighbors[k];
                    if(!Phydro[index]->Active) continue;
#ifdef __SKIP_HOTGAS__
                    if(Phydro[index]->UPred*Pall.ConvertUtoT >= 1.e+4) continue;
#endif
                    Phydro[index]->DuCooling += (1.e+4*Pall.ConvertTtoU-Phydro[index]->U)/Phydro[index]->dt_hydro;
                    Phydro[index]->U = Phydro[index]->UPred = Pall.ConvertTtoU*1.e+4;
#ifdef __EXPORT_TIMESTEP__
                    //Phydro[index]->k_hydro_localmin = MIN(Phydro[index]->k_hydro_localmin,PstarBody(i)->k-__DT_DIFF__);
                    Phydro[index]->k_hydro_localmin = MIN(Phydro[index]->k_hydro_localmin,PstarBody(i)->k+__DT_DIFF__);
#endif // __EXPORT_TIMESTEP__
#ifdef TASK_TEST_STROMGRENSPHERE
                    HIIFlag[index] = true;
#endif //TASK_TEST_STROMGRENSPHERE
                }
            }while(CurrentNodeID != RootNodeID);
        }
    }

    MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
    MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

    for(int i=0;i<NImport;i++){ // for imported data
        int nlist = 0;
        int RootNodeID = 0;
        int CurrentNodeID = HydroNode[RootNodeID].Children;
        do {
            CurrentNodeID = GetNeighborsIterativeApproach(CurrentNodeID,
                    HIIExportRecv[i].Pos,HIIExportRecv[i].Radius,&nlist,Neighbors);
            for(int k=0;k<nlist;k++){
                int index = Neighbors[k];
                if(!Phydro[index]->Active) continue;
#ifdef __SKIP_HOTGAS__
                if(Phydro[index]->UPred*Pall.ConvertUtoT >= 1.e+4) continue;
#endif
                Phydro[index]->DuCooling += (10000.0*Pall.ConvertTtoU-Phydro[index]->U)/Phydro[index]->dt_hydro;
                Phydro[index]->U = Phydro[index]->UPred = Pall.ConvertTtoU*1.e+4;
#ifdef __EXPORT_TIMESTEP__
                //Phydro[index]->k_hydro_localmin = MIN(Phydro[index]->k_hydro_localmin,HIIExportRecv[i].k_star);
                Phydro[index]->k_hydro_localmin = MIN(Phydro[index]->k_hydro_localmin,HIIExportRecv[i].k_star+__DT_DIFF__);
#endif // __EXPORT_TIMESTEP__
#ifdef TASK_TEST_STROMGRENSPHERE
                HIIFlag[index] = true;
#endif //TASK_TEST_STROMGRENSPHERE
            }
        }while(CurrentNodeID != RootNodeID);
    }

#ifdef __EXPORT_TIMESTEP__
    double EraMinimum = Pall.EraLocal+0.1*Pall.dtnow;
    for(int i=0;i<Pall.Nhydro;i++){
        if(Phydro[i]->Active == false){
            //if(Phydro[i]->k_hydro - Phydro[i]->k_hydro_localmin > __DT_DIFF__){
                //double dt_localmin = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin+__DT_DIFF__));
            if(Phydro[i]->k_hydro > Phydro[i]->k_hydro_localmin){
                double dt_localmin = Pall.dtmin*exp2((double)(Phydro[i]->k_hydro_localmin));
                int step = (int)(Pall.EraLocal/dt_localmin);
                double NextUpdateEra;
                do{
                    NextUpdateEra = step*dt_localmin;
                    step ++;
                }while(NextUpdateEra < EraMinimum);
                Phydro[i]->NextUpdateEra = NextUpdateEra;
            }
        }
    }
#endif // __EXPORT_TIMESTEP__

    return ;
}


/*
 * This fuction calculates the sizes of HII regions around ``active'' young
 * stars. The evaluation method for these sizes is similar to that for the
 * kernel size determination of SPH particles. 
 */
static void CalcHIIregions(void){

    double t = GetElapsedTime();

    int MyID = MPIGetMyID();
    int NProcs = MPIGetNumProcs();
    int Neighbors[MaxNeighborSize];
    MPI_Status  mpi_status;

    int NActiveYoungStars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->HIIflag){
#ifdef PRESERVE_SNII_EVENTRATE // {
            double MassInMsun = Pstar[i]->InitialMass*Pall.UnitMass/MSUN_CGS;
            double prob = SNIINumber*MassInMsun;
            if(prob >= 1.0){
                NActiveYoungStars ++;
            } else {
                if(Pstar[i]->TypeIIProb){
                    NActiveYoungStars ++;
                } else {
                    Pstar[i]->HIIflag = false;
                }
            }
#else //PRESERVE_SNII_EVENTRATE //} //{
            NActiveYoungStars ++;
#endif // PRESERVE_SNII_EVENTRATE //}
        }
    }
    // dprintlmpi(NActiveYoungStars);

    static int HIIExportFlagsMaxAllocated = 0;

    static bool *HIIExportFlags;
    static struct StructActiveHIIParticle *ActiveHIIParticle;

    if(HIIExportFlagsMaxAllocated < MAX(NActiveYoungStars,NAdditionUnit)){
        if(HIIExportFlagsMaxAllocated > 0){
            free(HIIExportFlags);
            free(ActiveHIIParticle);
        }
        HIIExportFlagsMaxAllocated = (int)(MAX(ForAngelsShare*NActiveYoungStars,NAdditionUnit));

        HIIExportFlags = malloc(sizeof(bool)*HIIExportFlagsMaxAllocated*NProcs);
        ActiveHIIParticle = malloc(sizeof(struct StructActiveHIIParticle)*HIIExportFlagsMaxAllocated);
    }

    int log = NActiveYoungStars;
    NActiveYoungStars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->HIIflag){
            int Offset = NActiveYoungStars*NProcs;
            HIIExportFlags[Offset+NProcs-1] = true;
            ActiveHIIParticle[NActiveYoungStars].Index = i;
            ActiveHIIParticle[NActiveYoungStars].Nlist = 0;

            // Set the initial Radius.
            if(Pstar[i]->StromgrenRadius>0.e0){ // Reuse old one.
                ActiveHIIParticle[NActiveYoungStars].Radius = Pstar[i]->StromgrenRadius;
            } else{ // Set an initial guess.
                ActiveHIIParticle[NActiveYoungStars].Radius = 0.25*PstarBody(i)->Eps;
            }
#ifdef __PHOTON_COUNT_BASE__
            ActiveHIIParticle[NActiveYoungStars].PhotonCount = 0.e0;
#else //__PHOTON_COUNT_BASE__
            ActiveHIIParticle[NActiveYoungStars].Mass = 0.e0;
#endif // __PHOTON_COUNT_BASE__

#ifdef PRESERVE_SNII_EVENTRATE
            double MassInMsun = Pstar[i]->InitialMass*Pall.UnitMass/MSUN_CGS;
            double prob = SNIINumber*MassInMsun;
            if(prob >= 1.0){
#   ifdef USE_ASRFLX_NLY // {
                double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[i]->FormationTime)/YEAR_CGS;
                ActiveHIIParticle[NActiveYoungStars].LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                    ASRFLXGetNLY(Age,Pstar[i]->Z);
#   else //USE_ASRFLX_NLY //}//{
                ActiveHIIParticle[NActiveYoungStars].LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                    ReturnNumberofLyPhoton(i,0);
#   endif // USE_ASRFLX_NLY // }
            } else {
                if(Pstar[i]->TypeIIProb){
#   ifdef USE_ASRFLX_NLY // {
                    double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[i]->FormationTime)/YEAR_CGS;
                    ActiveHIIParticle[NActiveYoungStars].LyAphoton = (1.0/SNIINumber)*
                        ASRFLXGetNLY(Age,Pstar[i]->Z);
#   else //USE_ASRFLX_NLY //}//{
                    ActiveHIIParticle[NActiveYoungStars].LyAphoton = (1.0/SNIINumber)*ReturnNumberofLyPhoton(i,0);
#   endif // USE_ASRFLX_NLY // }
                } 
            }
#else
#   ifdef USE_ASRFLX_NLY // {
            double Age = Pall.UnitTime*(Pall.TCurrent - Pstar[i]->FormationTime)/YEAR_CGS;
            ActiveHIIParticle[NActiveYoungStars].LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                ASRFLXGetNLY(Age,Pstar[i]->Z);
#   else //USE_ASRFLX_NLY //}//{
            ActiveHIIParticle[NActiveYoungStars].LyAphoton = (PstarBody(i)->Mass*Pall.UnitMass/MSUN_CGS)*
                ReturnNumberofLyPhoton(i,0);
#   endif // USE_ASRFLX_NLY // }
#endif
            ActiveHIIParticle[NActiveYoungStars].Pos[0] = PstarBody(i)->PosP[0];
            ActiveHIIParticle[NActiveYoungStars].Pos[1] = PstarBody(i)->PosP[1];
            ActiveHIIParticle[NActiveYoungStars].Pos[2] = PstarBody(i)->PosP[2];
            ActiveHIIParticle[NActiveYoungStars].Rvalue = ActiveHIIParticle[NActiveYoungStars].Lvalue = 0.e0;

            ActiveHIIParticle[NActiveYoungStars].HIIRegion = true;

            NActiveYoungStars ++;
        }
    }
    assert(log == NActiveYoungStars);


    int BitMask = 0x01; 
    int NExportThisTime[NProcs-1];
    int NImportThisTime[NProcs-1];
    int NExportThisTimeNew[NProcs];
    int NImportThisTimeNew[NProcs];

    struct StructHIIExport *HIIExportSend[NProcs-1];
    struct StructHIIExport *HIIExportRecv = NULL;
    struct StructHIIImport *HIIImportSend = NULL;
    struct StructHIIImport *HIIImportRecv[NProcs-1];
    MPI_Status mpi_status_Export_Send[NProcs-1];
    MPI_Request mpi_request_Export_Send[NProcs-1];
    MPI_Status mpi_status_Export_Recv[NProcs-1];
    MPI_Request mpi_request_Export_Recv[NProcs-1];

    int GlobalNActiveYoungStars;
    MPI_Allreduce(&NActiveYoungStars,&GlobalNActiveYoungStars,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    // if(MyID == MPI_ROOT_RANK)
        // dprintlmpi(GlobalNActiveYoungStars);

    CheckLocalDomainEdge();

    // Node center
    double BoxCenter[] = {GravityNode[0].Pos[0],GravityNode[0].Pos[1],GravityNode[0].Pos[2]};

    int Niteration = 0;
    do{
        for(int i=0;i<NActiveYoungStars;i++){ 
            int Offset = i*NProcs;
            if(HIIExportFlags[Offset+NProcs-1]){
                ActiveHIIParticle[i].DistMin =  0.e0;
#ifdef __PHOTON_COUNT_BASE__
                ActiveHIIParticle[i].PhotonCount = 0.e0;
                ActiveHIIParticle[i].PhotonCountDistMin = 0.e0;
#else //__PHOTON_COUNT_BASE__
                ActiveHIIParticle[i].Mass = 0.e0;
                ActiveHIIParticle[i].MassDistMin = 0.e0;
#endif // __PHOTON_COUNT_BASE__
                ActiveHIIParticle[i].Nlist = 0;
                for(int k=0;k<NProcs-1;k++)
                    HIIExportFlags[Offset+k] = false;
                LocalKernelMax = fmax(LocalKernelMax,
                        DISTANCE(BoxCenter,ActiveHIIParticle[i].Pos)
                            +2.0*ActiveHIIParticle[i].Radius);
            }
        }

        int NExportMaxThisTime = 0;
        for(int i=0;i<NProcs-1;i++){ // Check Export
            NExportThisTime[i] = CheckHIIExportFlagsModified(i,NProcs,HIIExportFlags,NActiveYoungStars,ActiveHIIParticle);
            CheckSizeofBufferExportSendIndex(NExportThisTime[i],sizeof(struct StructHIIExport),i);
            CheckSizeofBufferImportRecvIndex(NExportThisTime[i],sizeof(struct StructHIIImport),i);
            HIIExportSend[i] = BufferExportSend[i];
            HIIImportRecv[i] = BufferImportRecv[i];

            int NExport = 0;
            if(NExportThisTime[i]>0){
            for(int k=0;k<NActiveYoungStars;k++){
                int Offset = k*NProcs;
                if(HIIExportFlags[Offset+NProcs-1]){ 
                    if(HIIExportFlags[Offset+i]&BitMask){ 
                        HIIExportSend[i][NExport].Pos[0] = ActiveHIIParticle[k].Pos[0];
                        HIIExportSend[i][NExport].Pos[1] = ActiveHIIParticle[k].Pos[1];
                        HIIExportSend[i][NExport].Pos[2] = ActiveHIIParticle[k].Pos[2];
                        HIIExportSend[i][NExport].Radius = ActiveHIIParticle[k].Radius;
                        HIIExportSend[i][NExport].Leaf = k;
                        NExport ++;
                    }
                }
            }
            }
            NExportThisTime[i] = NExport;
        }

        int NImport = 0;
        int NImportThisTime2[NProcs];
        int NExportThisTime2[NProcs];
        NImportThisTime2[MPIGetMyID()] = 0;
        for(int i=0;i<NProcs-1;i++){
            NExportThisTime2[CommunicationTable[i].SendRank] = NExportThisTime[i];
        }
        MPI_Alltoall(NExportThisTime2,1,MPI_INT,NImportThisTime2,1,MPI_INT,MPI_COMM_WORLD);
        for(int i=0;i<NProcs-1;i++){
            NImportThisTime[i] = NImportThisTime2[CommunicationTable[i].RecvRank];
            NImport += NImportThisTime[i];
        }
        int NImportAll = NImport;

        CheckSizeofBufferExportRecv(NImportAll,sizeof(struct StructHIIExport));
        CheckSizeofBufferImportSend(NImportAll,sizeof(struct StructHIIImport));
        HIIExportRecv = BufferExportRecv;
        HIIImportSend = BufferImportSend; 

        NImport = 0;
        int counter_send = 0;
        int counter_recv = 0;
        int SendFlag,RecvFlag;
        for(int i=0;i<NProcs-1;i++){
            if(NExportThisTime[i]>0){
                MPI_Isend(HIIExportSend[i],
                    NExportThisTime[i]*sizeof(struct StructHIIExport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_HII_DENSITY_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NImportThisTime[i]>0){
                MPI_Irecv(HIIExportRecv+NImport,
                    NImportThisTime[i]*sizeof(struct StructHIIExport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_HII_DENSITY_EXPORT+i,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
                NImport += NImportThisTime[i];
            }
        }

        for(int i=0;i<NActiveYoungStars;i++){  // Check local
            int Offset = i*NProcs;
            if(HIIExportFlags[Offset+NProcs-1]){
                struct StructHIILocalInfo TemporalData =
                    ReturnNeighborsInfo(ActiveHIIParticle+i,Neighbors);

                ActiveHIIParticle[i].Nlist = TemporalData.Nlist;
#ifdef __PHOTON_COUNT_BASE__
                ActiveHIIParticle[i].PhotonCount = TemporalData.PhotonCount;
                ActiveHIIParticle[i].PhotonCountDistMin = TemporalData.PhotonCountDistMin;
#else //__PHOTON_COUNT_BASE__
                ActiveHIIParticle[i].Mass = TemporalData.Mass;
                ActiveHIIParticle[i].MassDistMin = TemporalData.MassDistMin;
#endif // __PHOTON_COUNT_BASE__
                ActiveHIIParticle[i].DistMin = TemporalData.DistMin;

#ifdef __LOCAL__
                ActiveHIIParticle[i].LocalUpdateFlag = false;
                int IsLocal = 0;
                for(int k=0;k<NProcs-1;k++){
                    if(HIIExportFlags[Offset+k]&BitMask){
                        IsLocal ++;
                    }
                }
                if(IsLocal == 0){
                    ActiveHIIParticle[i].LocalUpdateFlag = true;
                    UpdateHIIRadiusLocalModified(i,NActiveYoungStars,ActiveHIIParticle,Neighbors,
                            MyID,NProcs,HIIExportFlags);
                }
#endif //__LOCAL__
            }
        }

        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);


        for(int i=0;i<NImportAll;i++){
#if 0
            struct StructActiveHIIParticle TempActiveHIIParticle = {
                .Pos[0]=HIIExportRecv[i].Pos[0],
                .Pos[1]=HIIExportRecv[i].Pos[1],
                .Pos[2]=HIIExportRecv[i].Pos[2],
                .Radius=HIIExportRecv[i].Radius,
                .DistMin=0,
            };
            struct StructHIILocalInfo TemporalData =
                ReturnNeighborsInfo(&TempActiveHIIParticle,Neighbors);

#endif
#if 1
            struct StructHIILocalInfo TemporalData =
                ReturnNeighborsInfo(&(struct StructActiveHIIParticle)
                    {.Pos[0]=HIIExportRecv[i].Pos[0],
                     .Pos[1]=HIIExportRecv[i].Pos[1],
                     .Pos[2]=HIIExportRecv[i].Pos[2],
                     .Radius=HIIExportRecv[i].Radius,
                     .DistMin=0,},Neighbors);
#endif

            HIIImportSend[i].Nlist = TemporalData.Nlist;
            HIIImportSend[i].DistMin = TemporalData.DistMin;
#ifdef __PHOTON_COUNT_BASE__
            HIIImportSend[i].PhotonCount = TemporalData.PhotonCount;
            HIIImportSend[i].PhotonCountDistMin = TemporalData.PhotonCountDistMin;
#else //__PHOTON_COUNT_BASE__
            HIIImportSend[i].Mass = TemporalData.Mass;
            HIIImportSend[i].MassDistMin = TemporalData.MassDistMin;
#endif // __PHOTON_COUNT_BASE__
            HIIImportSend[i].Leaf = HIIExportRecv[i].Leaf;
        }

        NImportAll = 0;
        int NImportAllNew = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew[i] = 0;
            for(int k=0;k<NImportThisTime[i];k++){
                if(HIIImportSend[NImportAll].Nlist > 0){
                    HIIImportSend[NImportAllNew] = HIIImportSend[NImportAll];
                    NImportThisTimeNew[i] ++;
                    NImportAllNew ++;
                }
                NImportAll ++;
            }
        }

        int NImportThisTimeNew2[NProcs];
        int NExportThisTimeNew2[NProcs];
        NImportThisTimeNew2[MPIGetMyID()] = 0;
        for(int i=0;i<NProcs-1;i++){
            NImportThisTimeNew2[CommunicationTable[i].SendRank] = NImportThisTimeNew[i];
        }
        MPI_Alltoall(NImportThisTimeNew2,1,MPI_INT,NExportThisTimeNew2,1,MPI_INT,MPI_COMM_WORLD);
        for(int i=0;i<NProcs-1;i++){
            NExportThisTimeNew[i] = NExportThisTimeNew2[CommunicationTable[i].RecvRank];
        }


        NImport = 0;
        counter_send = counter_recv = 0;
        for(int i=0;i<NProcs-1;i++){
            if(NImportThisTimeNew[i]>0){
                MPI_Isend(HIIImportSend+NImport,
                    NImportThisTimeNew[i]*sizeof(struct StructHIIImport),
                        MPI_BYTE,CommunicationTable[i].SendRank,TAG_HII_DENSITY_IMPORT,
                            MPI_COMM_WORLD,mpi_request_Export_Send+counter_send);
                MPI_Test(mpi_request_Export_Send+counter_send,&SendFlag,MPI_STATUS_IGNORE);
                counter_send ++;
            }
            if(NExportThisTimeNew[i]>0){
                MPI_Irecv(HIIImportRecv[i],
                    NExportThisTimeNew[i]*sizeof(struct StructHIIImport),
                        MPI_BYTE,CommunicationTable[i].RecvRank,TAG_HII_DENSITY_IMPORT,
                            MPI_COMM_WORLD,mpi_request_Export_Recv+counter_recv);
                MPI_Test(mpi_request_Export_Recv+counter_recv,&RecvFlag,MPI_STATUS_IGNORE);
                counter_recv ++;
            }
            NImport += NImportThisTimeNew[i];
        }
        MPI_Waitall(counter_send,mpi_request_Export_Send,mpi_status_Export_Send);
        MPI_Waitall(counter_recv,mpi_request_Export_Recv,mpi_status_Export_Recv);

        for(int i=0;i<NProcs-1;i++){
            for(int k=0;k<NExportThisTimeNew[i];k++){ 
                int leaf = HIIImportRecv[i][k].Leaf;
                ActiveHIIParticle[leaf].Nlist += HIIImportRecv[i][k].Nlist;
#ifdef __PHOTON_COUNT_BASE__
                ActiveHIIParticle[leaf].PhotonCount += HIIImportRecv[i][k].PhotonCount;
#else // __PHOTON_COUNT_BASE__
                ActiveHIIParticle[leaf].Mass += HIIImportRecv[i][k].Mass;
#endif //__PHOTON_COUNT_BASE__
                if((ActiveHIIParticle[leaf].DistMin > HIIImportRecv[i][k].DistMin)&&(HIIImportRecv[i][k].Nlist > 0)){
                    ActiveHIIParticle[leaf].DistMin = HIIImportRecv[i][k].DistMin;
#ifdef __PHOTON_COUNT_BASE__
                    ActiveHIIParticle[leaf].PhotonCountDistMin = HIIImportRecv[i][k].PhotonCountDistMin;
#else // __PHOTON_COUNT_BASE__
                    ActiveHIIParticle[leaf].MassDistMin = HIIImportRecv[i][k].MassDistMin;
#endif //__PHOTON_COUNT_BASE__
                }
            }
        }

        int LocalNActiveYoungStars = CheckLocalMassAndUpdateHIIRadiusModified2(NActiveYoungStars,
                NProcs,HIIExportFlags,ActiveHIIParticle);

        MPI_Allreduce(&LocalNActiveYoungStars,&GlobalNActiveYoungStars,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

#if 1
        //if(MyID == MPI_ROOT_RANK) dprintlmpi(GlobalNActiveYoungStars); fflush(NULL); MPI_Barrier(MPI_COMM_WORLD);

        //if(GlobalNActiveYoungStars < 10){
        if(NActiveYoungStars < 10){
            int NLocalActiveLeaves = 0;
            for(int i=0;i<NActiveYoungStars;i++){
                int Offset = i*NProcs;
                // if(HIIExportFlags[Offset+NProcs-1]){ 
#ifdef __PHOTON_COUNT_BASE__ //{
                    double Qest = ActiveHIIParticle[i].PhotonCount;
#else // __PHOTON_COUNT_BASE__ //{ //}
                    double Qest = fact*(SQ(ActiveHIIParticle[i].Mass*Pall.UnitMass)/(CUBE(2.0*ActiveHIIParticle[i].Radius*Pall.UnitLength)));
#endif // __PHOTON_COUNT_BASE__ //}

                    double Qmin = (1.0-ConvergenceFactor)*ActiveHIIParticle[i].LyAphoton;
                    double Qmax = (1.0+ConvergenceFactor)*ActiveHIIParticle[i].LyAphoton;

                // }
            }
        }
#endif

        Niteration ++;
        if(Niteration > 100)
            break;
    } while (0<GlobalNActiveYoungStars);

    if(MPIGetMyID()==MPI_ROOT_RANK){
        fprintf(stderr,"%d iterations for Stromgren radius determination.\n",Niteration);
        fprintf(stderr,"SR update time is %g [sec]\n",GetElapsedTime()-t);
    }

    for(int i=0;i<NActiveYoungStars;i++){
        int Index = ActiveHIIParticle[i].Index;
#ifdef __PHOTON_COUNT_BASE__
        Pstar[Index]->Density = ActiveHIIParticle[i].PhotonCount;
#else // __PHOTON_COUNT_BASE__
        Pstar[Index]->Density = ActiveHIIParticle[i].Mass;
#endif // __PHOTON_COUNT_BASE__
        if(ActiveHIIParticle[i].HIIRegion){
            Pstar[Index]->StromgrenRadius = 2*ActiveHIIParticle[i].Radius;
        } else {
            Pstar[Index]->StromgrenRadius = NONE;
        }
    }

    return ;
}

static bool HIIfirst = true;

void InitializeHIIregions(void){

    if(HIIfirst == true){
        LoadnLyData();
        AllocateContanctedDomainID();
        HIIfirst = false;
    }
    return ;
}

void HIIregions(void){

#ifndef USE_HIIREGION_MODEL
    return ;
#endif // USE_HIIREGION_MODEL
    if(Pall.NActivesStars_t == 0)
        return ;
    // MPI_Barrier(MPI_COMM_WORLD);

    if(HIIfirst == true){
        LoadnLyData();
        AllocateContanctedDomainID();
        HIIfirst = false;
    }

#ifdef EVALUATE_SIZES_ALL_TOGETHER //{

    int NLocalActiveYoungStars = ReturnCalcSizeElementNumber(CS_TypeHII,false);
    if(ReturnCalcSizeElementNumber(CS_TypeHII,true) == 0){
        MakeHIIregion();
        return ;
    }
#else // EVALUATE_SIZES_ALL_TOGETHER //}//{

    // Check HII flags.
    int NLocalActiveYoungStars = CheckHIIflags(0);
    int GlobalActiveYoungStars;

    MPI_Allreduce(&NLocalActiveYoungStars,&GlobalActiveYoungStars,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if(GlobalActiveYoungStars == 0){
        MakeHIIregion();
        return ;
    }

    CheckContactedDomain();
    CalcHIIregions();
#endif // EVALUATE_SIZES_ALL_TOGETHER //}

    // Write HII raidus.
    // WriteHIIradius();

    //UpdateLocalFUVvalue();

    int NActiveYoungStars = 0, GlobalNActiveYoungStars = 0;
    for(int i=0;i<Pall.Nstars;i++){
        if(Pstar[i]->HIIflag)
            NActiveYoungStars ++;
    }
    MPI_Allreduce(&NActiveYoungStars,&GlobalNActiveYoungStars,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(GlobalNActiveYoungStars> 0){
        MakeHIIregion();
    }

    return ;
}

