/* StellarWind.c */
void UpdateStellarWindEnergy(const int index);
void InitStellarWind(void);
double GetStellarWindSourceLifeTime(void);
double GetStellarWindMaxTimestep(void);
bool CheckStellarWindSource(const int index);
bool CheckStellarWindSourceThroughPbody(const int index);
int CountStellarWindSourceNumber(void);
