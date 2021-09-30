/* TestRuns.c */
#ifdef TASK_TEST_GRAVITYTREE
int main_GravityTreeTest(const int argc, char *argv[]);
#endif //TASK_TEST_GRAVITYTREE
#ifdef TASK_TEST_HYDROTREE
int main_HydroTreeTest(const int argc, char *argv[]);
int main_HydroTreeIntegralTest(const int argc, char *argv[]);
int main_HydroTreeRobustTest(const int argc, char *argv[]);
#endif //TASK_TEST_HYDROTREE
#ifdef TASK_TEST_NEIGHBORSEARCH
int main_Test_NeighborSearch(const int argc, char *argv[]);
#endif //TASK_TEST_NEIGHBORSEARCH
#ifdef TASK_TEST_HYDRO_QUANTITIES //{
int main_Test_HydroQuantities(const int argc, char *argv[]);
#endif // TASK_TEST_HYDRO_QUANTITIES //}

#ifdef TASK_TEST_SINKPARTICLE
int main_Test_Sink(const int argc, char *argv[]);
int main_Test_SinkParticleRun(const int argc, char *argv[]);
#endif //TASK_TEST_SINKPARTICLE

#ifdef TASK_TEST_FOF //{
int main_FOFTest(const int argc, char *argv[]);
#endif // TASK_TEST_FOF //{

#ifdef TASK_TEST_STELLARFEEDBACK // TASK_TEST_STELLARFEEDBACK //{
int main_StellarFeedbackTest(const int argc, char *argv[]);
#endif // TASK_TEST_STELLARFEEDBACK //}

#ifdef TASK_TEST_MOMENTUMFEEDBACK // TASK_TEST_MOMENTUMFEEDBACK //{
int main_MomentumFeedbackTest(const int argc, char *argv[]);
#endif // TASK_TEST_MOMENTUMFEEDBACK //}

#ifdef TASK_TEST_EQUILIBRIUM_TEMPERATURE //{
int main_Test_EquilibriumTemperature(const int argc, char *argv[]);
#endif // TASK_TEST_EQUILIBRIUM_TEMPERATURE //}

#ifdef TASK_TEST_1D_THERMAL_CONDUCTIVITY //{
int main_Test_1DThermalConductivity(const int argc, char *argv[]);
#endif // TASK_TEST_1D_THERMAL_CONDUCTIVITY //}

#ifdef TASK_TEST_3D_THERMAL_CONDUCTIVITY //{
int main_Test_3DThermalConductivity(const int argc, char *argv[]);
#endif // TASK_TEST_3D_THERMAL_CONDUCTIVITY //}

#ifdef TASK_TEST_FUVFEEDBACK //{
int main_TestFUVFeedback(const int argc, char *argv[]);
#endif // TASK_TEST_FUVFEEDBACK //}

#ifdef TASK_TEST_STROMGRENSPHERE
int main_TestStromgrenSphere(const int argc, char *argv[]);
#endif //TASK_TEST_STROMGRENSPHERE

#ifdef TASK_TEST_RADIATION_PRESSURE //{
int main_TestRadiationPressure(const int argc, char *argv[]);
#endif // TASK_TEST_RADIATION_PRESSURE //}

#ifdef TASK_TEST_DOMAINDECOMPOSITION //{
int main_TestDomainDecomposition(const int argc, char *argv[]);
#endif // TASK_TEST_DOMAINDECOMPOSITION //}

#ifdef TASK_TEST_PARTICLE_SPLITTING //{
int main_Test_ParticleSplitting(const int argc, char *argv[]);
#endif // TASK_TEST_PARTICLE_SPLITTING //}
