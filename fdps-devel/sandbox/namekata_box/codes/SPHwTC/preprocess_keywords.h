#define Current_Code_Mode (1)
#define Normal_Operation_Mode (0)
#define Make_Glass_Mode (1)
#define Shock_Tube_Test_Mode (2)
#define Blastwave_Test_Mode (3)
#define Surface_Tension_Test_Mode (4)
#define KHI_Test_Mode (5)
#define RTI_Test_Mode (6)
#define Blob_Test_Mode (7)
#define Gresho_Vortex_Test_Mode (8)
#define Noh_Problem_Mode (9)
#define Keplerian_Ring_Test_Mode (10)
#define Evrard_Test_Mode (11)
#define Restart_Mode (0)
#define Make_Rate_Tables (0)
#define Simulation_Dimension (3)
#define Coupled_with_Hydro (1)
#define Consider_Size_Distribution (0)
#define Multicomponent_ISM (1)
#define Bound_Gas_Temperature (1)
#define Use_Temperature_Floor (1)
#define Temperature_Floor_Type (2)
#define Floor_By_EpsG (1)
#define Floor_By_MSPH (2)
#define Treatment_EnergyEq (2)
#define Adiabatic (1)
#define Isothermal (2)
#define Barotropic (3)
#define Cooling (4)
#define Radtr_Cooling (5)
#define Radiation_Force (1)
#define Update_Chemical_Abundance (1)
#define Update_Gas_Temperature (1)
#define Update_Dust_Temperature (1)
#define Thermal_Conduction (0)
#define Thermal_Conductivity_Model (1)
#define Constant_Conductivity (1)
#define Spitzer_Conductivity (2)
#define Global_Timestepping (1)
#define Timescale_Dyn (1)
#define Timescale_Engy (2)
#define Timescale_Chem (3)
#define External_Gravity (0)
#define Self_Gravity (1)
#define Force_Error_Analysis (0)
#define SPH_Formulation (3)
#define Standard_SPH (1)
#define Springel2005 (2)
#define Hopkins2013_Energy (3)
#define Hopkins2013_Entropy (4)
#define Saitoh_Makino_2012 (5)
#define SPH_Kernel (2)
#define M4_Cubic_Spline (1)
#define M4_CS_w_TC92 (2)
#define M5_Quartic_Spline (3)
#define M6_Quintic_Spline (4)
#define Core_Triangle (5)
#define HOCT4 (6)
#define Wendland_C2 (7)
#define Wendland_C4 (8)
#define Wendland_C6 (9)
#define Artificial_Viscosity (2)
#define Monaghan_Gingold_1983 (1)
#define Monaghan_1997 (2)
#define Balsara_Switch (0)
#define Morris_Monaghan_Switch (0)
#define Cullen_Dehnen_Switch (0)


#pragma once
//=============================
// Conditionally-defined macros
//=============================
#if (Current_Code_Mode == Make_Glass_Mode)
#ifdef Simulation_Dimension
#undef Simulation_Dimension
#endif
#define Simulation_Dimension (3)
#endif

#if (Current_Code_Mode == Make_Glass_Mode)
#ifdef Coupled_with_Hydro
#undef Coupled_with_Hydro
#endif
#define Coupled_with_Hydro (1)
#endif

#if (Current_Code_Mode == Make_Glass_Mode)
#ifdef Treatment_EnergyEq
#undef Treatment_EnergyEq
#endif
#define Treatment_EnergyEq (2)
#endif

#if (Current_Code_Mode == Make_Glass_Mode)
#ifdef External_Gravity
#undef External_Gravity
#endif
#define External_Gravity (0)
#endif

#if (Current_Code_Mode == Make_Glass_Mode)
#ifdef Self_Gravity
#undef Self_Gravity
#endif
#define Self_Gravity (0)
#endif

#if (Current_Code_Mode == Shock_Tube_Test_Mode)
#ifdef Simulation_Dimension
#undef Simulation_Dimension
#endif
#define Simulation_Dimension (3)
#endif

#if (Current_Code_Mode == Shock_Tube_Test_Mode)
#ifdef Multicomponent_ISM
#undef Multicomponent_ISM
#endif
#define Multicomponent_ISM (0)
#endif

#if (Current_Code_Mode == Shock_Tube_Test_Mode)
#ifdef Bound_Gas_Temperature
#undef Bound_Gas_Temperature
#endif
#define Bound_Gas_Temperature (0)
#endif

#if (Current_Code_Mode == Shock_Tube_Test_Mode)
#ifdef Treatment_EnergyEq
#undef Treatment_EnergyEq
#endif
#define Treatment_EnergyEq (1)
#endif

#if (Current_Code_Mode == Shock_Tube_Test_Mode)
#ifdef Global_Timestepping
#undef Global_Timestepping
#endif
#define Global_Timestepping (1)
#endif

#if (Current_Code_Mode == Shock_Tube_Test_Mode)
#ifdef External_Gravity
#undef External_Gravity
#endif
#define External_Gravity (0)
#endif

#if (Current_Code_Mode == Shock_Tube_Test_Mode)
#ifdef Self_Gravity
#undef Self_Gravity
#endif
#define Self_Gravity (0)
#endif



//=============================
// Error conditions            
//=============================
#if ((Simulation_Dimension != 1) && (Simulation_Dimension != 2) && (Simulation_Dimension != 3))
#error
#endif

#if (((Treatment_EnergyEq == Cooling) || (Treatment_EnergyEq == Radtr_Cooling)) && ((Update_Chemical_Abundance == 0) && (Update_Gas_Temperature == 0) && (Update_Dust_Temperature == 0)))
#error
#endif

#if (Update_Chemical_Abundance == 0)
#error
#endif

#if ((Multicomponent_ISM == 0) && ((Treatment_EnergyEq == Cooling) ||  (Treatment_EnergyEq == Radtr_Cooling)))
#error
#endif

#if ((Balsara_Switch == 1) && (Cullen_Dehnen_Switch == 1))
#error
#endif

#if ((Morris_Monaghan_Switch == 1) && (Cullen_Dehnen_Switch == 1))
#error
#endif

