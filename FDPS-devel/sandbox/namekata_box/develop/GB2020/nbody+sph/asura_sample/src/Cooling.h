/* Cooling.c */
void InitializeCoolingTable(void);
void CalcCooling(void);
double MetalInterpolatedCoolingRate(const double T, const double metal);

// For Spaans (2008)
void InitializeCoolingRateLookUpTable(void);
void CalcCoolingSpaans2008(void);
double InterpolatedCoolingRateSpaans2008(double Tg, double Zg, double G0g, double fH2g);
