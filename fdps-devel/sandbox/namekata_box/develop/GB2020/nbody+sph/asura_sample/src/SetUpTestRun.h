/* SetUpTestRun.c */
void ActivateAllparticles(void);
void SetViscosityParameters(const double AlphaMin, const double AlphaMax, const double SourceFactor, const double InvDumpingScale);

void Init3DCollapseTest(const int number);
void InitSelfSimilarCooling(const int number);
void ExternalPotentialForSelfSimilarCooling(void);
void Init3DShockTubeGrid(const int number);
void Init3DShockTube(const int number);
void Read3DShockTube(const int number, char LoadFileName[]);
void GlassCondition3DShockTube(void);
void OutPut3DCollapseTest(char *name);
void FileIO3DCollapse(void);
void FileIOBlastWave(void);
void OutPutZeroPlane(char *name);
void FileIOSelfSimilar(void);
void InitSphereWithMeanParticleDistance(const int Number);
double CalcTotalEnergyForColdCollapse(void);
void OutPutColdCollapse(char *fname);
void InitializeNavarroWhiteTest(const int NGrid);
void OutPutDarkMatter(char *fname_dm);
void OutPutNavarroWhite(char *fname_hydro, char *fname_star, char *fname_dm);
void InitTreeTest(const int Number);
void InitIsothermalSphericalCollapseTest(const int number);
void InitStandardIsothermalTestCase(const int number);
void FileOutputIsothermalSphericalCollapse(void);
void InitBreakingDam(const int MultiplicationConstant);
void GravitationalForceFromEarth(void);
void AntiForceFromBounaryRegions(const int MultiplicationConstant);
void InitUniformSphereTest(const int number);
void InitRandomSphereTest(const int number);
void InitUniformBoxTest(const int number);
void InitRandomBoxTestHydro(const int number);
void InitRandomBoxTest(const int number);
void ReadKeisokuCheckParallel(const int number);
void ReadKeisokuBenchParallel(void);
// void ReadMultiMassCosmologicalInitialCondition(char fbody[], char fboundary[]);
// void ReadMultiMassNbodyInitialCondition(void);
void OutPutNFW(void);
void ParticleDistributionGenetator(int number, unsigned long int seed);
void OutPutAllParticlesInASCIIFormat(void);

#ifdef TASK_MERGER
void ReadGalactICS(void);
void ReadBinaryGalactICS(char fname[], const double GasFlactionDisk, const double GasFlactionHalo);
void ReadBinaryGalactICSFlag(char fname[]);
void ReadBinaryMergingGalaxies(char fname[]);
#endif //TASK_MERGER
#ifdef TASK_DICE_RUN //{
void ReadDiceASCII(char fname[]);
void ReadDiceBinary(char fname[]);
#endif // TASK_DICE_RUN  //}
#ifdef TASK_MW
void InitializeParticleDistributionForMilkyWay(const int Number, const double MassFactor);
void InitializeParticleDistributionForMilkyWayWithInvRProf(const int Number, const double MassFactor);
void InitializeParticleDistributionForMilkyWayWithHaloGas(const int Number);
void MilkyWayPotentialHaloDisk(void);
double CircularVelocityMilkyWayHaloDisk(const double r);
void ParticleRemover(const double Radius);
void CheckKennicutLaw(const int mode);
#endif //TASK_MW
#ifdef TASK_CLOUD_EXPLOSION
void ReadInitCloudExplosion(void);
void InitCloudExplosion(const int number);
void AddThermalEnergy(const double fraction);
#endif //TASK_CLOUD_EXPLOSION
#ifdef TASK_BLAST_WAVE
void InitBlastWave(const int number);
void ReadBlastWave(void);
void InitBlastWaveRandom(const int number);
void AddSmoothedEnergyBlastWave(const int number);
void AddSmoothedEnergyBlastWaveDirectSearch(const int PeakIndex, const int number);
void GlassCondition(const double RhoTrue);
void OutPutBlastWaveAll(char *name);
#endif //TASK_BLAST_WAVE
#ifdef TASK_SINUSOIDAL_WAVE 
void InitSinusoidalWave(void);
#endif // TASK_SINUSOIDAL_WAVE 
#ifdef TASK_NFW
void ReadNFWInitialCondition(void);
void ReadNFWTdyn(void);
void ExternalPotentialForNFW(void);
#endif //TASK_NFW
#ifdef TASK_M2_COLLAPSE
void InitM2SphericalCollapseTest(const int number);
#endif //TASK_M2_COLLAPSE
#ifdef TASK_TURBULENCE
void InitTurbulentCloud(char fname[], const int mode, const double LeftMetal, const double RightMetal);
#endif //TASK_TURBULENCE
#ifdef TASK_ROTATINGDISK_WITH_SINK
void InitRotatingDiskWithSink(const int Number);
#endif //TASK_ROTATINGDISK_WITH_SINK
#ifdef TASK_AGNTORUS
void InitAGNTorus(const int Number, const int DstortionFlag);
void AGNForceFromExternalPotentials(void);
#endif //TASK_AGNTORUS
#ifdef TASK_GALACTIC_CENTER
void InitGCCloud(double InitPos[], double InitVel[]);
#endif // TASK_GALACTIC_CENTER
#ifdef TASK_SANTABARBARA
void ReadSantaBarbaraInitialCondition(const int mode);
#endif // TASK_SANTABARBARA
#ifdef TASK_GALAXY_FORMATION
void ReadMusic(char fname[], const bool BinaryFlag);
#endif // TASK_GALAXY_FORMATION


#ifdef TASK_COLD_COLLAPSE
void InitColdCollapseTest(const int Number, const double rv, const int MakePotential);
void InitColdCollapseTestMixed(const int Number, const double rv, const double fraction, const int MakePotential);
#endif //TASK_COLD_COLLAPSE
#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
void InitSphericalParticleDistributionWithTwoTypes(const int Number, const double fraction, const double Radius, const double Distance);
void RestoreSphericalParticleDistribution(const int Number, double Pos[][3], double Mass[], double Eps[], const double M_p, const double Eps_p, const double Radius, const double Distance);
void InitSphericalParticleDistributionForErrorEstimation(const int Number, const double fraction, const double Radius);
void RestoreSphericalParticleDistributionErrorEstimation(const int Number, double Pos[][3], double Mass[], double Eps[]);
#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR
#ifdef TASK_TEST_STROMGRENSPHERE
void InitTestStromgrenSphere(const double nH, const double Radius, const int Nhydro, const int Nstars, double Pos[restrict][3], double LyAlpha[restrict]);
#endif //TASK_TEST_STROMGRENSPHERE
#ifdef TASK_1D_SHOCKE_TUBE
void InitShockTube(const int Number);
void InitSodShockTube(const int Number);
void InitShockTubeHK(const int Number);
void InitShockTube123Problem(const int Number);
void OutPutShockTube(char *fname);
#endif //TASK_1D_SHOCKE_TUBE
#ifdef TASK_HYDROSTATIC
void InitHydroStatic(const int NGrid, const int FlagEqualMass, const double DensityRatio);
#endif // HYDROSTATIC
#ifdef TASK_KELVINHELMHOLZ_INSTABILITY
void InitKH(const int NGrid, const int mode, const int seed, const double density_contrast, const double TEnd);
#endif // TASK_KELVINHELMHOLZ_INSTABILITY
#ifdef TASK_TEST_1D_THERMAL_CONDUCTIVITY //{
void Init1DThermalConductivity(const int Number, const double T_JUMP);
#endif // TASK_TEST_1D_THERMAL_CONDUCTIVITY //}
#ifdef TASK_TEST_3D_THERMAL_CONDUCTIVITY //{
void Init3DThermalConductivity(const int Number, const int mode, const double T_JUMP);
#endif // TASK_TEST_3D_THERMAL_CONDUCTIVITY //}
#ifdef TASK_1D_TWOFLUIDS //{
void Init1DTwoFluids(const int Number);
#endif // TASK_1D_TWOFLUIDS //}
#ifdef TASK_KEPLER //{
void InitKepler(const int Number);
void KeplerPotential(void);
#endif // TASK_KEPLER //}
