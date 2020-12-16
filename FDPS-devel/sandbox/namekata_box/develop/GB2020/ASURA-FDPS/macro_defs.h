#pragma once

//#####################################
// Macro constants for SPH calculation
//#####################################
// (1) Type of SPH kernel (macro ASURA_FDPS_SPH_KERNEL)
#define ASURA_FDPS_BSPLINE_M4  (0)
#define ASURA_FDPS_BSPLINE_M5  (1)
#define ASURA_FDPS_BSPLINE_M6  (2)
#define ASURA_FDPS_WENDLAND_C2 (3)
#define ASURA_FDPS_WENDLAND_C4 (4)
#define ASURA_FDPS_WENDLAND_C6 (5)
// Check the correctness of the correspoinding compile-time macro
#if !defined(ASURA_FDPS_SPH_KERNEL)
#error Macro ASURA_FDPS_SPH_KERNEL is not defined.
#endif
#if ((ASURA_FDPS_SPH_KERNEL != ASURA_FDPS_BSPLINE_M4) && \
     (ASURA_FDPS_SPH_KERNEL != ASURA_FDPS_BSPLINE_M5) && \
     (ASURA_FDPS_SPH_KERNEL != ASURA_FDPS_BSPLINE_M6) && \
     (ASURA_FDPS_SPH_KERNEL != ASURA_FDPS_WENDLAND_C2) && \
     (ASURA_FDPS_SPH_KERNEL != ASURA_FDPS_WENDLAND_C4) && \
     (ASURA_FDPS_SPH_KERNEL != ASURA_FDPS_WENDLAND_C6))
#error The value of macro ASURA_FDPS_SPH_KERNEL is wrong.
#endif
// (2) Options for kernel size determination
// (3) Options for prediction method in kernel size determination (macro ASURA_FDPS_KERNEL_SIZE_PRED_METHOD)
#define ASURA_FDPS_THACKER_2000          (0)
#define ASURA_FDPS_THACKER_2000_MODIFIED (1)
// Check the correctness of the corresponding compile-time macro
#if !defined(ASURA_FDPS_KERNEL_SIZE_PRED_METHOD)
#define ASURA_FDPS_KERNEL_SIZE_PRED_METHOD (1) // the default method
#else
#if ((ASURA_FDPS_KERNEL_SIZE_PRED_METHOD != ASURA_FDPS_THACKER_2000) && \
     (ASURA_FDPS_KERNEL_SIZE_PRED_METHOD != ASURA_FDPS_THACKER_2000_MODIFIED))
#error The value of macro ASURE_FDPS_KERNEL_SIZE_PRED_METHOD is wrong.
#endif
#endif
// (4) Options for the evaluation of grad-h term (macro ASURA_FDPS_GRAD_H_TERM)
#define ASURA_FDPS_EQ15_IN_HOPKINS_2013 (0)
#define ASURA_FDPS_EQ18_IN_HOPKINS_2013 (1)
// Check the correctness of the corresponding compile-time macro
#if !defined(ASURA_FDPS_GRAD_H_TERM)
#define ASURA_FDPS_GRAD_H_TERM (1) // the default grad-h term
#else
#if ((ASURA_FDPS_GRAD_H_TERM != ASURA_FDPS_EQ15_IN_HOPKINS_2013) && \
     (ASURA_FDPS_GRAD_H_TERM != ASURA_FDPS_EQ18_IN_HOPKINS_2013))
#error The value of macro ASURA_FDPS_GRAD_H_TERM is wrong.
#endif
#endif
// (5) Options to control the artificial viscosity (macro ASURA_FDPS_ARTIFICIAL_VISCOSITY)
#define ASURA_FDPS_CONSTANT_VISCOSITY   (0)
#define ASURA_FDPS_MORRIS_MONAGHAN_1997 (1)
#define ASURA_FDPS_CULLEN_DEHNEN_2010   (2) // NOT YET AVAILABLE
// Check the correctness of the corresponding compile-time macro
#if !defined(ASURA_FDPS_ARTIFICIAL_VISCOSITY)
#define ASURA_FDPS_ARTIFICIAL_VISCOSITY (1) // the default artificial viscosity model
#else
#if ((ASURA_FDPS_ARTIFICIAL_VISCOSITY != ASURA_FDPS_CONSTANT_VISCOSITY) && \
     (ASURA_FDPS_ARTIFICIAL_VISCOSITY != ASURA_FDPS_MORRIS_MONAGHAN_1997) && \
     (ASURA_FDPS_ARTIFICIAL_VISCOSITY != ASURA_FDPS_CULLEN_DEHNEN_2010))
#error The value of macro ASURA_FDPS_ARTIFICIAL_VISCOSITY is wrong.
#endif
#endif

//#########################################
// Macro constants for gravity calculation
//#########################################
 

//#############################################
// Macro constants for special execution modes 
//#############################################
#define ASURA_FDPS_NORMAL_RUN                 (0)
#define ASURA_FDPS_GLASS_DATA_GENERATION_MODE (1)
#define ASURA_FDPS_PARTICLE_COLLISION_TEST    (2)
#define ASURA_FDPS_SHOCK_TUBE_TEST            (3)
#define ASURA_FDPS_SURFACE_TENSION_TEST       (4)
#define ASURA_FDPS_KHI_TEST                   (5) // Kelvin-Helmholtz instability test; NOT AVAILABLE YET
#define ASURA_FDPS_RTI_TEST                   (6) // Rayleigh-Taylor instability test; NOT AVAILABLE YET
#define ASURA_FDPS_POINT_EXPLOSION_TEST       (7) 
// Check the correctness of the correspoinding compile-time macro
#if !defined(ASURA_FDPS_EXEC_MODE)
#error Macro ASURA_FDPS_EXEC_MODE is not defined.
#endif
#if ((ASURA_FDPS_EXEC_MODE != ASURA_FDPS_NORMAL_RUN) && \
     (ASURA_FDPS_EXEC_MODE != ASURA_FDPS_GLASS_DATA_GENERATION_MODE) && \
     (ASURA_FDPS_EXEC_MODE != ASURA_FDPS_PARTICLE_COLLISION_TEST) && \
     (ASURA_FDPS_EXEC_MODE != ASURA_FDPS_SHOCK_TUBE_TEST) && \
     (ASURA_FDPS_EXEC_MODE != ASURA_FDPS_SURFACE_TENSION_TEST) && \
     (ASURA_FDPS_EXEC_MODE != ASURA_FDPS_KHI_TEST) && \
     (ASURA_FDPS_EXEC_MODE != ASURA_FDPS_RTI_TEST) && \
     (ASURA_FDPS_EXEC_MODE != ASURA_FDPS_POINT_EXPLOSION_TEST))
#error The value of macro ASURA_FDPS_EXEC_MODE is wrong.
#endif

// Automatic setting of necessary macros
//----------------------------------
// (1) Glass data generation
#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_GLASS_DATA_GENERATION_MODE
// turn on the use of isothermal EOS
#ifndef ASURA_FDPS_USE_ISOTHERMAL_EOS
#define ASURA_FDPS_USE_ISOTHERMAL_EOS
#endif
// turn on velocity dumping
#ifndef ASURA_FDPS_ENABLE_VELOCITY_DUMPING
#define ASURA_FDPS_ENABLE_VELOCITY_DUMPING
#endif
// turn off gravity calculation
#ifdef ASURA_FDPS_ENABLE_GRAVITY
#undef ASURA_FDPS_ENABLE_GRAVITY
#endif
// turn off cooling/heating
#ifdef ASURA_FDPS_ENABLE_COOLING_HEATING
#undef ASURA_FDPS_ENABLE_COOLING_HEATING
#endif
// turn off star formation
#ifdef ASURA_FDPS_ENABLE_STAR_FORMATION
#undef ASURA_FDPS_ENABLE_STAR_FORMATION
#endif
// turn off stellar feedback
#ifdef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#undef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#endif
#endif // ASURA_FDPS_GLASS_DATA_GENERATION_MODE
//----------------------------------
// (2) Particle collision test
#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_PARTICLE_COLLISION_TEST
// enforce fixed kernel size
#ifndef ASURA_FDPS_USE_FIXED_KERNEL_SIZE
#define ASURA_FDPS_USE_FIXED_KERNEL_SIZE
#endif
// turn off the use of isothermal EOS
#ifdef ASURA_FDPS_USE_ISOTHERMAL_EOS
#undef ASURA_FDPS_USE_ISOTHERMAL_EOS
#endif
// turn off gravity calculation
#ifdef ASURA_FDPS_ENABLE_GRAVITY
#undef ASURA_FDPS_ENABLE_GRAVITY
#endif
// turn off cooling/heating
#ifdef ASURA_FDPS_ENABLE_COOLING_HEATING
#undef ASURA_FDPS_ENABLE_COOLING_HEATING
#endif
// turn off star formation
#ifdef ASURA_FDPS_ENABLE_STAR_FORMATION
#undef ASURA_FDPS_ENABLE_STAR_FORMATION
#endif
// turn off stellar feedback
#ifdef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#undef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#endif
#endif // ASURA_FDPS_PARTICLE_COLLISION_TEST
//--------------------------
// (3) Shock tube test
#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_SHOCK_TUBE_TEST
// turn off the use of isothermal EOS
#ifdef ASURA_FDPS_USE_ISOTHERMAL_EOS
#undef ASURA_FDPS_USE_ISOTHERMAL_EOS
#endif
// turn off gravity calculation
#ifdef ASURA_FDPS_ENABLE_GRAVITY
#undef ASURA_FDPS_ENABLE_GRAVITY
#endif
// turn off cooling/heating
#ifdef ASURA_FDPS_ENABLE_COOLING_HEATING
#undef ASURA_FDPS_ENABLE_COOLING_HEATING
#endif
// turn off star formation
#ifdef ASURA_FDPS_ENABLE_STAR_FORMATION
#undef ASURA_FDPS_ENABLE_STAR_FORMATION
#endif
// turn off stellar feedback
#ifdef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#undef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#endif
#endif // ASURA_FDPS_SHOCK_TUBE_TEST
//--------------------------
// (4) Surface Tension Test
#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_SURFACE_TENSION_TEST
// check if macro PARTICLE_SIMULATOR_TWO_DIMENSION is defined or not.
#if !defined(PARTICLE_SIMULATOR_TWO_DIMENSION)
#error Compile-time macro PARTICLE_SIMULATOR_TWO_DIMENSION is not defined.
#endif
// turn off the use of isothermal EOS
#ifdef ASURA_FDPS_USE_ISOTHERMAL_EOS
#undef ASURA_FDPS_USE_ISOTHERMAL_EOS
#endif
// turn off gravity calculation
#ifdef ASURA_FDPS_ENABLE_GRAVITY
#undef ASURA_FDPS_ENABLE_GRAVITY
#endif
// turn off cooling/heating
#ifdef ASURA_FDPS_ENABLE_COOLING_HEATING
#undef ASURA_FDPS_ENABLE_COOLING_HEATING
#endif
// turn off star formation
#ifdef ASURA_FDPS_ENABLE_STAR_FORMATION
#undef ASURA_FDPS_ENABLE_STAR_FORMATION
#endif
// turn off stellar feedback
#ifdef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#undef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#endif
#endif // ASURA_FDPS_SURFACE_TENSION_TEST
//--------------------------
// (5) KHI test
//--------------------------
// (6) RTI test
//--------------------------
// (7) Point explosion test
#if ASURA_FDPS_EXEC_MODE == ASURA_FDPS_POINT_EXPLOSION_TEST
// turn off the use of isothermal EOS
#ifdef ASURA_FDPS_USE_ISOTHERMAL_EOS
#undef ASURA_FDPS_USE_ISOTHERMAL_EOS
#endif
// turn off gravity calculation
#ifdef ASURA_FDPS_ENABLE_GRAVITY
#undef ASURA_FDPS_ENABLE_GRAVITY
#endif
// turn off cooling/heating
#ifdef ASURA_FDPS_ENABLE_COOLING_HEATING
#undef ASURA_FDPS_ENABLE_COOLING_HEATING
#endif
// turn off star formation
#ifdef ASURA_FDPS_ENABLE_STAR_FORMATION
#undef ASURA_FDPS_ENABLE_STAR_FORMATION
#endif
// turn off stellar feedback
#ifdef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#undef ASURA_FDPS_ENABLE_STELLAR_FEEDBACK
#endif
#endif // ASURA_FDPS_POINT_EXPLOSION_TEST
//--------------------------


//##############################
// Mathematical macro functions
//##############################
#define SQ(x)   ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#define PWR2(x) ((x)*(x))
#define PWR3(x) ((x)*(x)*(x))
#define PWR4(x) ((x)*(x)*(x)*(x))
#define PWR5(x) ((x)*(x)*(x)*(x)*(x))
#define PWR6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define PWR7(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x))
#define PWR8(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))
#define PWR9(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))

#define SIGN(x) (((x)>0)-((x)<0))
