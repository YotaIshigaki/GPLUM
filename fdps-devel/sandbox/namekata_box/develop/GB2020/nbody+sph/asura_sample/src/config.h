#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <tgmath.h>
#include <assert.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#ifdef USE_CELIB //{
#include <CELib.h>
#endif // USE_CELIB //}
//#ifdef USE_FUVFEEDBACK //{
#include <ASRFLX.h>
//#endif // USE_FUVFEEDBACK //}
#ifdef USE_ASRCH //{
#include <ASRCH.h>
#endif // USE_ASRCH //}
#include "Astro.h"

//#define DBG_ID  (6818)
//#define DBG_ID  (146)
// #define DBG_ID  (3965)
//#define DBG_ID  (98825)
//#define DBG_ID  (61066)
//#define DBG_ID  (-1)
//#define DBG_ID  (95247)
//#define DBG_ID  (161196)
//#define DBG_ID  (450204)
//#define DBG_ID  (44445)
//#define DBG_ID  (95458)
//#define DBG_ID  (95980)
//#define DBG_ID  (4)
//#define DBG_ID  (888)
//#define DBG_ID  (284063)
#define DBG_ID  (2204569)

#define MaxIteration        (200)
#define DecompNSample   (1000) // The number of Sample particles broadcast for other nodes.
#define NAdditionUnit   (1000) // Minimum increase of particles at a moment.

#define ForAngelsShare          (1.2)   // For Angel's Share. Allocate reserved memory.
#define MinimumAllocationSize   (10)    // Minimum Allocation Size for malloc.
#define FirstAllocationSize    (100)

/***** G6A GRAPE-5 Compatible mode *****/
#if (!defined(HAVE_GRAPE5)&&!defined(HAVE_GRAPE7)&&!defined(HAVE_PHANTOM_GRAPE)&&!defined(HAVE_AVX_PHANTOM_GRAPE))
#define JMEMSIZE    (1024*64-1)
#endif
/***** G6A GRAPE-5 Compatible mode *****/

/***** Run Type *****/
#ifdef COSMOLOGICAL_RUN
#define CosmologicalRun (ON)
#else
#define CosmologicalRun (OFF)
#endif
/***** Run Type *****/

/***** Parameters For Tree *****/
//#define TreeOpeningAngle        (0.5)
//#define TreeOpeningAngle        (0.5)
#define TreeOpeningAngle        (0.5)

#define TreeNodeGenerationLimitNumberForGrav    (32)
//#define TreeNodeGenerationLimitNumberForGrav    (1)
#define TreeNodeGenerationLimitNumberForNBS     (8)
#define TreeMaxNodeLevel	    (31) // When long long int (or long int) is used for Key.
#define TreeSofteningFactor     (2.0)

#if (defined(HAVE_PHANTOM_GRAPE)||defined(HAVE_AVX_PHANTOM_GRAPE))
#define TreeNGroup              (128)
//#define TreeNGroup              (8)
//#define TreeNGroup              (2)
#elif defined(HAVE_GRAPE7) 
#define TreeNGroup              (512)
//#define TreeNGroup              (1024)
//#define TreeNGroup              (2048)
#else
#define TreeNGroup              (512)
//#define TreeNGroup              (10000)
#endif
/***** Parameters For Tree *****/

/***** Parameters For Neighbor Search *****/
#define MaxNeighborSize         (MAX_NEIGHBOR_SIZE)
/***** Parameters For Neighbor Search *****/

/***** Parameters For Individual Time Step *****/
//#define VariableTimeStep    (ON)
#ifdef USE_VARIABLE_TIMESTEP
//#ifdef USE_INDIVIDUAL_TIMESTEP
//#else
//#define MaximumTimeHierarchy    (0)
//#endif
#else 
//#if (VariableTimeStep == OFF)
//#define DT_CONST    (0.005*GIGAYEAR_CGS/Pall.UnitTime)
//#define DT_CONST    (1.e-3)
//#define DT_CONST    (exp2(-9.0))
#define DT_CONST    (exp2(-3.0))
#endif
#define MaximumTimeHierarchy    (MAXIMUM_TIME_HIERARCHY)
/***** Parameters For Individual Time Step *****/

/***** For time integlation *****/
#define TFactorDiffAcc    (0.1)
//#if MaximumTimeHierarchy == 0
#ifndef MAXIMUM_TIME_HIERARCHY
#define TFactor    (0.3)
#define TFactorCourant  (0.3)
#define TFactorU        (0.3)
#define TFactorVel      (0.5)
#define TFactorAcc      (0.1)
#else
//#define TFactor    (0.25)
#define TFactor    (0.3)
#define TFactorCourant  (0.3)
#define TFactorU        (0.3)
#define TFactorVel      (0.5)
#define TFactorAcc      (0.1)
/*
#define TFactor    (0.05)
#define TFactorCourant  (0.125)
#define TFactorAcc      (0.05)
*/
#endif
#define TFactorDecomposition    (0.01)
/***** For time integlation *****/


/***** Parameters For Hydro *****/
#define UseCooling          (ON)
#define SingleStepCoolingSkip   (SKIP_COOLING)
#define HeliumAbandance   (0.24) // Helium abandance.
// HydroTimeStep Version
#define HydroViscNormal (OFF)
#define HydroViscVelSig (ON)
/***** Parameters For Hydro *****/


/***** Parameters and Models For Starformation *****/
#define UseSFModelConvert   (OFF)
#define UseSFModelSpawn     (ON)
#define MaxSpawnTimes       (3)

#define UseConstantSFrate   (OFF)

#define IMFTYPE_SP          (1)

#define SFeff               (SF_EFFICIENCY)
#define SFDensityCrit       (SFCONDITION_DENSITY_CRITERION)
#define SFTemperatureCrit   (SFCONDITION_TEMPERATURE_CRITERION)
/***** Parameters and Models For Starformation *****/


/***** Types of Particles *****/
enum {
    TypeDM,
    TypeHydro,
    TypeStar,
    TypeSink,
    NTypes,
};
/***** Types of Particles *****/

#define ENUMTOSTR(var) #var


#define DIMENSION           (3)
#define BisectionDimension  (3)
#define TINY                (1.e-30)

/***** For Run Status *****/
enum {
    NewSimulation,
    RestartSimulation,
    DataFileOperation_Merge,
    DataFileOperation_Halve,
    DataFileOperation_Double,
};
/***** For Run Status *****/

/***** For Header Flags *****/
enum{
    HD_SizeofHeader,
    HD_SizeofInt,
    HD_CompactFormatFlag,
    HD_CompactDoubleFormatFlag,
    HD_CompactMixFormatFlag,
    HD_LeanFormatFlag,
};
/***** For Header Flags *****/

/***** For I/O *****/
//#define DATA_DUMP_INTERVAL   (30.0) // in unit of minute.
//#define ASCIIDATA_DUMP_INTERVAL   (10.0) // in unit of minute.
/***** For I/O *****/

#include "reconfig.h"
#include "Astro.h"
#include "Constants.h"
#include "Unit.h"
#include "Errors.h"
#include "debug.h"
#include "IO.h"
#include "DataStructures.h"
#include "BufferOperation.h"
#include "StructureOperation.h"
#include "ParallelOperation.h"
#include "PeriodicWrapping.h"
#include "Allocators.h"
#include "MPIParameters.h"
#include "MPITags.h"
#include "Utilities.h"
#include "ElapsedTime.h"
#include "Cosmological.h"
#include "CosmologicalTime.h"
#include "ParticleAddRemoveOperation.h"
#include "RandomNumberGenerator.h"
#include "GRAPEEmulator.h"
#include "CheckStructures.h"
#include "BoundaryCondition.h"
#include "PowderSnow.h"

