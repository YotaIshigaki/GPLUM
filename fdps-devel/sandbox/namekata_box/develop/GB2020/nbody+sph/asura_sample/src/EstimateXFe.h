#pragma once 

enum {
    XFeFeH_SAGAFull,    // SAGA database, MW+dwarfs
    XFeFeH_SAGACompact, // SAGA database, MW
    XFeFeH_GALAH,       // GALAH DR2.1 data, flag_cannon = 0 and flag_xfe = 0
    XFeFeH_APOGEE,      // APOGEE data, snr>500
    XFeFeH_GALAH_APOGEE, // Some of missing primary elements in GALAH data is covered by APOGEE data.
    XFeFeH_GALAH_SAGAFull, // Some of missing primary elements in GALAH data is covered by SAGA full data.
    XFeFeH_Number,
};


double EstimateXFe(const char ElementName[], const double FeH, const int DataType);
double EstimateNX(const char ElementName[], const double FeH, const int DataType);
double EstimateMassX(const char ElementName[], const double FeH, const int DataType);

void EstimateAbudancePattern(double Elements[], const double FeH, const int DataType);
void EstimateAbudancePatternForMass(double Elements[], const double FeH, const double Mass, const int DataType);
void EstimateAbundancePatternTest(const int DataType);
