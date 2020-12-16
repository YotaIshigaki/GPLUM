#include <assert.h>

#if DEBUG_LEVEL >= 10
#define ASSERT(x)   assert(x)
#else
#define ASSERT(x)   /* NOP */
#endif

#if DEBUG_LEVEL >= 30
#define LOG(args...)   fprintf(stderr,args)
#else
#define LOG(args...)   /* NOP */
#endif

#define dbg(...) \
    (printf("dbg: %s %u @%s:",__FILE__,__LINE__,__func__), \
    printf(" "__VA_ARGS__)) 

#define TargetGlobalIDBody(x) \
    for(int i=0;i<Pall.Ntotal;i++) if(Pbody[i]->GlobalID == x)

#define TargetGlobalIDHydro(x) \
    for(int i=0;i<Pall.Nhydro;i++) if(PhydroBody(i)->GlobalID == x)

#define TargetGlobalIDStar(x) \
    for(int i=0;i<Pall.Nstars;i++) if(PstarBody(i)->GlobalID == x)
