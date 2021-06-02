#include "config.h"

/** Element Pointers **/
struct  StructPall     Pall;
struct  StructTimingResults     TimingResults;
StructPbodyptr    PbodyElements = NULL;   // Linked List
StructPhydroptr   PhydroElements = NULL;  // Linked List
StructPstarptr    PstarElements = NULL;   // Linked List
StructPsinkptr    PsinkElements = NULL;   // Linked List
unsigned int PbodyElementsSize = 0;
unsigned int PhydroElementsSize = 0;
unsigned int PstarElementsSize = 0;
unsigned int PsinkElementsSize = 0;
/** Element Pointers **/

/** Direct Access Pointer **/
StructPbodyptr *Pbody;
StructPhydroptr *Phydro;
StructPstarptr *Pstar;
StructPsinkptr *Psink;
unsigned int PbodySize = 0;
unsigned int PhydroSize = 0;
unsigned int PstarSize = 0;
unsigned int PsinkSize = 0;
/** Direct Access Pointer **/

/* Communication Buffer */
int *NumberofBufferExportSendAllocated;
int NumberofBufferExportRecvAllocated;
int NumberofBufferImportSendAllocated;
int *NumberofBufferImportRecvAllocated;
void **BufferExportSend;
void *BufferExportRecv;
void *BufferImportSend;
void **BufferImportRecv;
/* Communication Buffer */

/* Hydro Interaction Buffer */
int NumberofBufferHydroInteractionFlagsAllocated;
bool *BufferHydroInteractionFlags;
/* Hydro INteraction  Buffer */

/* Communication Buffer */
struct  StructCommunicationTable   *CommunicationTable;
/* Communication Buffer */

/* Structures for Domain Decomposition */
struct StructInfoBiSection *InfoBiSection;
/* Structures for Domain Decomposition */

/* Edges of particle distributions */
struct StructEdges *EdgesForHydro = NULL;
struct StructEdges *EdgesForActiveHydro = NULL;
struct StructEdges *EdgesForGravity = NULL;
struct StructEdges *EdgesForStars = NULL;
struct StructEdges *EdgesForYoungStars = NULL;
struct StructEdges *EdgesForSink = NULL;
struct StructEdges *EdgesForActiveSink = NULL;
/* Edges of particle distributions */

/* Structures for Tree */
StructGravityRoot GravityRoot;
struct StructGravityNode *GravityNode;
struct StructGravityCache *GravityCache;
struct StructGravityAccPotCache *GravityAccPotCache;

StructHydroRoot HydroRoot;
struct StructHydroNode *HydroNode;
StructNBCache *NBCache;

StructHydroRoot HydroRootImported;
struct StructHydroNode *HydroNodeImported;
StructNBCache *NBImportedCache;

StructHydroRoot StellarRoot;
struct StructHydroNode *StellarNode;
StructNBCache *StellarNBCache;

StructHydroRoot SinkRoot;
struct StructHydroNode *SinkNode;
StructNBCache *SinkNBCache;

StructYoungStarRoot YoungStarRoot;
struct StructYoungStarNode *YoungStarNode;
struct StructYoungStarCache *YoungStarCache;
struct StructYoungStarResultCache *YoungStarResultCache;
/* Structures for Tree */


/* Friend-of-Friend */
int NFOFGroups;
int FOFSize = 0;
int FOFCatalogSize = 0;
struct StructFOF *FOF = NULL;
struct StructFOFCatalog *FOFCatalog = NULL;
/* Friend-of-Friend */

/* Log file pointers */
FILE *FpEnergy;
FILE *FpMomentum;
FILE *FpAMomentum;
FILE *FpElapsedTime;
FILE *FpStepLog;
char FnameEnergy[MaxCharactersInLine];
char FnameMomentum[MaxCharactersInLine];
char FnameAMomentum[MaxCharactersInLine];
char FnameElapsedTime[MaxCharactersInLine];
char FnameStepLog[MaxCharactersInLine];
/* Log file pointers */

/* Random Generator */
gsl_rng *RandomGenerator;
/* Random Generator */

/* MPI Communicators */
#ifdef USE_BARYON_COMM //{
MPI_Group hydro_group;
MPI_Group activehydro_group;
MPI_Comm MPI_HYDRO_COMM_WORLD = MPI_COMM_NULL;
MPI_Comm MPI_ACTIVEHYDRO_COMM_WORLD = MPI_COMM_NULL;
MPI_Group baryon_group;
MPI_Comm MPI_BARYON_COMM_WORLD = MPI_COMM_NULL;
MPI_Group activebaryon_group;
MPI_Comm MPI_ACTIVEBARYON_COMM_WORLD = MPI_COMM_NULL;
#endif // USE_BARYON_COMM //}
/* MPI Communicators */
