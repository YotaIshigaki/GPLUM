#include "config.h"

int GetNeighbors(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if(distance2 < hh){
                    list[nlist] = NBCache[leaf].Leaf;
                    nlist ++;
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}
	return nlist;
}

int GetNeighborsLimited(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if(distance2 < hh){
                    list[nlist] = NBCache[leaf].Leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}
	return nlist;
}

int GetActiveNeighborsLimited(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if((NBCache[leaf].Active == true)&&(distance2 < hh)){
                    list[nlist] = NBCache[leaf].Leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}
	return nlist;
}

int GetNeighborsIterativeApproach(const int StartNodeID, 
        double Pos[restrict], const double h, int *nlist, int list[restrict]){

    double h2 = SQ(h);
	*nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = StartNodeID; 
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
            if((*nlist)+HydroNode[CurrentNodeID].NumberofLeaves > MaxNeighborSize)
                return CurrentNodeID;
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if(distance2 < h2){
                    list[*nlist] = NBCache[leaf].Leaf;
                    (*nlist) ++;
                }
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Next; 
		}
	}
    return CurrentNodeID;
}

#ifdef USE_SMOOTHED_NEIGHBOR_NUMBER //{
static inline double KernelHydroKernel(const double r, const double InvKerneli) __attribute__((always_inline));
static inline double KernelHydroKernel(const double r, const double InvKerneli){ 

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

int GetNeighborsSmoothedNumberIterativeApproach(const int StartNodeID, 
        double Pos[restrict], const double h, int *nlist, int list[restrict], double *SmoothedNumber){

    double h2 = SQ(h);
	//*nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = StartNodeID; 
    double InvKerneli = 2.0/h;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
            if((*nlist)+HydroNode[CurrentNodeID].NumberofLeaves > MaxNeighborSize)
                return CurrentNodeID;
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if(distance2 < h2){
                    *SmoothedNumber += KernelHydroKernel(sqrt(distance2),InvKerneli);
                    list[*nlist] = NBCache[leaf].Leaf;
                    (*nlist) ++;
                }
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Next; 
		}
	}
    return CurrentNodeID;
}
#endif // USE_SMOOTHED_NEIGHBOR_NUMBER //{

int GetNeighborsIterativeApproachFOF(const int  StartNodeID, 
        double Pos[restrict], const double LinkingLength, int *nlist, int list[restrict]){

    double LinkingLength2 = SQ(LinkingLength);
	*nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = StartNodeID; 
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(LinkingLength+HydroNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
            if(*nlist+HydroNode[CurrentNodeID].NumberofLeaves > MaxNeighborSize)
                return CurrentNodeID;

			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if(distance2 < LinkingLength2){
                    list[*nlist] = NBCache[leaf].Leaf;
                    (*nlist) ++;
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}

    return CurrentNodeID;
}


int GetNeighborNumbers(double Pos[restrict], const double h){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if(distance2 < hh){
                    nlist ++;
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}
	return nlist;
}


int GetNeighborsPairs(double Pos[restrict], const double h, int list[restrict]){

    double hh = h*h;
    int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( (dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax))&&(dx2>SQ(HydroNode[CurrentNodeID].KernelMax)) ){
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        } else if(HydroNode[CurrentNodeID].Children != NONE){
            CurrentNodeID = HydroNode[CurrentNodeID].Children;
        } else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if( (distance2<hh) || (distance2<4.0*SQ(NBCache[leaf].Kernel)) ){
                    list[nlist] = NBCache[leaf].Leaf;
                    nlist ++;
                }
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        }
    }
    return nlist;
}

int GetNeighborsPairsLimited(double Pos[restrict], const double h, int list[restrict]){

    double hh = h*h;
    int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( (dx2>SQ(h+HydroNode[CurrentNodeID].DistanceMax))&&(dx2>SQ(HydroNode[CurrentNodeID].KernelMax)) ){
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        } else if(HydroNode[CurrentNodeID].Children != NONE){
            CurrentNodeID = HydroNode[CurrentNodeID].Children;
        } else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if( (distance2<hh) || (distance2<4.0*SQ(NBCache[leaf].Kernel)) ){
                    list[nlist] = NBCache[leaf].Leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        }
    }
    return nlist;
}

int GetNeighborsLimitedImported(double Pos[restrict], const double h, int list[restrict]){

    double hh = h*h;
    int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNodeImported[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNodeImported[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNodeImported[CurrentNodeID].Pos,Pos);
#endif
        if( dx2>SQ(h+HydroNodeImported[CurrentNodeID].DistanceMax)){
            CurrentNodeID = HydroNodeImported[CurrentNodeID].Next;
        } else if(HydroNodeImported[CurrentNodeID].Children != NONE){
            CurrentNodeID = HydroNodeImported[CurrentNodeID].Children;
        } else {
			int NumberofLeaves = HydroNodeImported[CurrentNodeID].NumberofLeaves;
            int header = HydroNodeImported[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                //if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBImportedCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBImportedCache[leaf].Pos,Pos);
#endif
                if(distance2 < hh){
                    list[nlist] = NBImportedCache[leaf].Leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
            }
            CurrentNodeID = HydroNodeImported[CurrentNodeID].Next;
        }
    }
    return nlist;
}

int GetNeighborsPairsLimitedImported(double Pos[restrict], const double h, int list[restrict]){

    double hh = h*h;
    int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNodeImported[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNodeImported[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNodeImported[CurrentNodeID].Pos,Pos);
#endif
        if( (dx2>SQ(h+HydroNodeImported[CurrentNodeID].DistanceMax))&&
                (dx2>SQ(HydroNodeImported[CurrentNodeID].KernelMax)) ){
            CurrentNodeID = HydroNodeImported[CurrentNodeID].Next;
        } else if(HydroNodeImported[CurrentNodeID].Children != NONE){
            CurrentNodeID = HydroNodeImported[CurrentNodeID].Children;
        } else {
			int NumberofLeaves = HydroNodeImported[CurrentNodeID].NumberofLeaves;
            int header = HydroNodeImported[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                //if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBImportedCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBImportedCache[leaf].Pos,Pos);
#endif
                if( (distance2<hh) || (distance2<4.0*SQ(NBImportedCache[leaf].Kernel)) ){
                    list[nlist] = NBImportedCache[leaf].Leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
            }
            CurrentNodeID = HydroNodeImported[CurrentNodeID].Next;
        }
    }
    return nlist;
}

int GetNeighborsPairsInteractiveApproach(const int StartNodeID, 
        double Pos[restrict], const double h, int *nlist, int list[restrict]){

    double hh = h*h;
    *nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = StartNodeID; 
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( (dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax))&&(dx2>SQ(HydroNode[CurrentNodeID].KernelMax)) ){
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        } else if(HydroNode[CurrentNodeID].Children != NONE){
            CurrentNodeID = HydroNode[CurrentNodeID].Children;
        } else {
            if((*nlist)+HydroNode[CurrentNodeID].NumberofLeaves > MaxNeighborSize)
                return CurrentNodeID;

			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if( (distance2<hh) || (distance2<4.0*SQ(NBCache[leaf].Kernel)) ){
                    list[*nlist] = NBCache[leaf].Leaf;
                    (*nlist) ++;
                }
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        }
    }
    return CurrentNodeID;
}

int GetNeighborsDirect(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;
	for(int i=0;i<Pall.Nhydro;i++){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],PhydroPosP(i)[k],k));
        }
#else
		double dx2 = DISTANCE2(Pos,PhydroPosP(i));
#endif
		if(dx2 <= hh){
			list[nlist] = i;
			nlist ++;
		}
	}

	return nlist;
}

int GetNeighborsPairsDirect(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;
	for(int i=0;i<Pall.Nhydro;i++){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],PhydroPosP(i)[k],k));
        }
#else
        double dx2 = DISTANCE2(Pos,PhydroPosP(i));
#endif
        if( (dx2<hh) || (dx2<4.0*SQ(Phydro[i]->Kernel)) ){
            list[nlist] = i;
            nlist ++;
		}
	}
	return nlist;
}

int GetNeighborsDirectLimitedImported(double Pos[restrict], const double h, int list[restrict]){

    double hh = h*h;
    int nlist = 0;

    int RootNodeID = 0;
    int NumberofLeaves = HydroNodeImported[RootNodeID].NumberofLeaves;
    int header = HydroNodeImported[RootNodeID].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
#ifdef PERIODIC_RUN 
        double distance2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            distance2 += SQ(PeriodicDistance(Pos[k],NBImportedCache[leaf].Pos[k],k));
        }
#else
        double distance2 = DISTANCE2(NBImportedCache[leaf].Pos,Pos);
#endif
        if( (distance2<hh) || (distance2<4.0*SQ(NBImportedCache[leaf].Kernel)) ){
            list[nlist] = NBImportedCache[leaf].Leaf;
            nlist ++;
            if(nlist==MaxNeighborSize)
                return nlist;
        }
    }
    return nlist;
}

int GetNeighborsPairsDirectLimitedImported(double Pos[restrict], const double h, int list[restrict]){

    double hh = h*h;
    int nlist = 0;

    int RootNodeID = 0;
    int NumberofLeaves = HydroNodeImported[RootNodeID].NumberofLeaves;
    int header = HydroNodeImported[RootNodeID].Leaves;
    for(int k=0;k<NumberofLeaves;k++){
        int leaf = header+k;
#ifdef PERIODIC_RUN 
        double distance2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            distance2 += SQ(PeriodicDistance(Pos[k],NBImportedCache[leaf].Pos[k],k));
        }
#else
        double distance2 = DISTANCE2(NBImportedCache[leaf].Pos,Pos);
#endif
        if( (distance2<hh) || (distance2<4.0*SQ(NBImportedCache[leaf].Kernel)) ){
            list[nlist] = NBImportedCache[leaf].Leaf;
            nlist ++;
            if(nlist==MaxNeighborSize)
                return nlist;
        }
    }
    return nlist;
}

int ReturnNeighborNumber(const int CurrentNodeID, double Pos[restrict], const double h){

    int Nlist = 0;
#ifdef PERIODIC_RUN 
    double dx2 = 0.e0;
    for(int k=0;k<DIMENSION;k++){
        dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
    }
#else
    double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif

    if(dx2 < SQ(h+HydroNode[CurrentNodeID].DistanceMax)){
		if(HydroNode[CurrentNodeID].Children == NONE){
	        double hh = h*h;
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++){
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
                }
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if(distance2 < hh){
                    Nlist ++;
                }
			}
        } else {
            int ChildNodeID = HydroNode[CurrentNodeID].Children;
            do{
                Nlist += ReturnNeighborNumber(ChildNodeID,Pos,h);
                ChildNodeID = HydroNode[ChildNodeID].Sister;
            } while (ChildNodeID != NONE);
        }
    }
    return Nlist;
}

int GetNeighborsLimitedNBCacheIndex(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+HydroNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		} else if(HydroNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = HydroNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if(distance2 < hh){
                    list[nlist] = leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
			}
			CurrentNodeID = HydroNode[CurrentNodeID].Next;
		}
	}
	return nlist;
}


int GetNeighborsPairsLimitedNBCacheIndex(double Pos[restrict], const double h, int list[restrict]){

    double hh = h*h;
    int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = HydroNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],HydroNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(HydroNode[CurrentNodeID].Pos,Pos);
#endif
        if( (dx2>SQ(h+HydroNode[CurrentNodeID].DistanceMax))&&(dx2>SQ(HydroNode[CurrentNodeID].KernelMax)) ){
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        } else if(HydroNode[CurrentNodeID].Children != NONE){
            CurrentNodeID = HydroNode[CurrentNodeID].Children;
        } else {
			int NumberofLeaves = HydroNode[CurrentNodeID].NumberofLeaves;
            int header = HydroNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(HydroRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],NBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(NBCache[leaf].Pos,Pos);
#endif
                if( (distance2<hh) || (distance2<4.0*SQ(NBCache[leaf].Kernel)) ){
                    list[nlist] = leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
            }
            CurrentNodeID = HydroNode[CurrentNodeID].Next;
        }
    }
    return nlist;
}

//////////////////////////////////////////////////////////

int GetNeighborsFromStellarTree(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = StellarNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],StellarNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(StellarNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+StellarNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = StellarNode[CurrentNodeID].Next;
		} else if(StellarNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = StellarNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = StellarNode[CurrentNodeID].NumberofLeaves;
            int header = StellarNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(StellarRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],StellarNBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(StellarNBCache[leaf].Pos,Pos);
#endif
                if(distance2 < hh){
                    list[nlist] = StellarNBCache[leaf].Leaf;
                    nlist ++;
                }
			}
			CurrentNodeID = StellarNode[CurrentNodeID].Next;
		}
	}
	return nlist;
}

int GetNeighborsLimitedFromStellarTree(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = StellarNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],StellarNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(StellarNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+StellarNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = StellarNode[CurrentNodeID].Next;
		} else if(StellarNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = StellarNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = StellarNode[CurrentNodeID].NumberofLeaves;
            int header = StellarNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(StellarRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],StellarNBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(StellarNBCache[leaf].Pos,Pos);
#endif
                if(distance2 < hh){
                    list[nlist] = StellarNBCache[leaf].Leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
			}
			CurrentNodeID = StellarNode[CurrentNodeID].Next;
		}
	}
	return nlist;
}

int GetActiveNeighborsLimitedFromStellarTree(double Pos[restrict], const double h, int list[restrict]){

	double hh = h*h;
	int nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = StellarNode[RootNodeID].Children;
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],StellarNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(StellarNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+StellarNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = StellarNode[CurrentNodeID].Next;
		} else if(StellarNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = StellarNode[CurrentNodeID].Children;
		} else {
			int NumberofLeaves = StellarNode[CurrentNodeID].NumberofLeaves;
            int header = StellarNode[CurrentNodeID].Leaves;
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(StellarRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],StellarNBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(StellarNBCache[leaf].Pos,Pos);
#endif
                if((StellarNBCache[leaf].Active == true)&&(distance2 < hh)){
                    list[nlist] = StellarNBCache[leaf].Leaf;
                    nlist ++;
                    if(nlist==MaxNeighborSize)
                        return nlist;
                }
			}
			CurrentNodeID = StellarNode[CurrentNodeID].Next;
		}
	}
	return nlist;
}

int GetNeighborsIterativeApproachFromStellarTree(const int StartNodeID, 
        double Pos[restrict], const double h, int *nlist, int list[restrict]){

    double h2 = SQ(h);
	*nlist = 0;

    int RootNodeID = 0;
    int CurrentNodeID = StartNodeID; 
	while(CurrentNodeID != RootNodeID){
#ifdef PERIODIC_RUN 
        double dx2 = 0.e0;
        for(int k=0;k<DIMENSION;k++){
            dx2 += SQ(PeriodicDistance(Pos[k],StellarNode[CurrentNodeID].Pos[k],k));
        }
#else
        double dx2 = DISTANCE2(StellarNode[CurrentNodeID].Pos,Pos);
#endif
        if( dx2 > SQ(h+StellarNode[CurrentNodeID].DistanceMax) ){
			CurrentNodeID = StellarNode[CurrentNodeID].Next;
		} else if(StellarNode[CurrentNodeID].Children != NONE){		
			CurrentNodeID = StellarNode[CurrentNodeID].Children;
		} else {
            if((*nlist)+StellarNode[CurrentNodeID].NumberofLeaves > MaxNeighborSize)
                return CurrentNodeID;
			int NumberofLeaves = StellarNode[CurrentNodeID].NumberofLeaves;
            int header = StellarNode[CurrentNodeID].Leaves;
            //dprintlmpi(NumberofLeaves);
            for(int k=0;k<NumberofLeaves;k++){
                int leaf = header+k;
                if(StellarRoot.Leaves[leaf] < 0) continue;
#ifdef PERIODIC_RUN 
                double distance2 = 0.e0;
                for(int l=0;l<DIMENSION;l++)
                    distance2 += SQ(PeriodicDistance(Pos[l],StellarNBCache[leaf].Pos[l],l));
#else
                double distance2 = DISTANCE2(StellarNBCache[leaf].Pos,Pos);
#endif
                if(distance2 < h2){
                    list[*nlist] = StellarNBCache[leaf].Leaf;
                    (*nlist) ++;
                }
            }
            CurrentNodeID = StellarNode[CurrentNodeID].Next; 
		}
	}
    return CurrentNodeID;
}
