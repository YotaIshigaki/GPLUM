#pragma once 

void InitializeDecomposition(void);
//void PreDomainDecomposition(void);
void PreDomainDecomposition(const int mode);
void InspectParticleFluctuation(void);
#if DECOMPOSITION_TYPE == 1 //{ 
int GetDomainKeyXYZ(double Pos[3]);
#endif // DECOMPOSITION_TYPE == 1 //}
