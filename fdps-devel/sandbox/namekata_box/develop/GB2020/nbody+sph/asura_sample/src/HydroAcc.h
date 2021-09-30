/* HydroAcc.c */
void CalcDuDtAcc(void);
void CalcDuDtAccForCorrection(void);
void CalcDuDtAccEnergyDensityForCorrection(void);
void CalcDuDtAccDirect(void);
void CalcSinkPressureForcesToHydro(double Pos[], double Vel[], double Mass, double Kernel, double Rho);
