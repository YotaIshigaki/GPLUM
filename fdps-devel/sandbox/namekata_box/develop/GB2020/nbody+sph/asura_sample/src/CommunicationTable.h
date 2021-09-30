#pragma once 

void InitializeCommunicationTable(void);

#ifdef USE_BARYON_COMM //{
void InitializeHydroCommunicationOrder(void);
void InitializeActiveHydroCommunicationOrder(void);
void InitializeBaryonCommunicationOrder(void);
void InitializeActiveBaryonCommunicationOrder(void);
#endif // USE_BARYON_COMM //}
