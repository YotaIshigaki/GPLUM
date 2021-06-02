#pragma once

double GetRadiationPressure(int index);
void UpdateBolometricLuminosity(const int index);
void InitRadiationPressure(void);
double GetRadiationPressureSourceLifeTime(void);
double GetRadiationPressureMaxTimestep(void);
bool CheckRadiationPressureSource(const int index);
bool CheckRadiationPressureSourceThroughPbody(const int index);
int CountRadiationPressureSourceNumber(void);
void RadiationPressureKick(void);
