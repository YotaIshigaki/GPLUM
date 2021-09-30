#!/usr/bin/perl -w
use strict;

# This script needs one argument from the command line.

# Check the number of the arguments
if ( @ARGV != 1 ) {
   die "The number of the arguments should be one.\n";
}
# Check the value of the second argument.
my $restart = $ARGV[0];
my $mktable;
if ( "$restart" eq "restart_run" ) {
   $restart = 1;
   $mktable = 0;
} elsif ( "$restart" eq "initial_run" ) {
   $restart = 0;
   $mktable = 0;
} elsif ( "$restart" eq "make_tables" ) {
   $restart = 0;
   $mktable = 1;
} else {
   die "The content of the second argument is invalid\n";
}

#-----------------------------------------------------------
# [1] Declare preprocess keywords and those values by hash
#-----------------------------------------------------------
our $pointer=-1;
our %pp_keywords;
our @ordering;

our $pointer_wc=-1;
our @pp_keywords_wc;
our @pp_values_wc;
our @pp_conditions_wc;
# `wc` is a short of `with conditions`

our $pointer_err=-1;
our @err_list;

#-----------------------------------------------------------
# [2] Define keywords and values
#-----------------------------------------------------------
# [!!!!! Important Notes !!!!!]
#   - Never change the keywords and values marked by #!.
#   - Macros with `future` are reserved macros for future implementation.
#
#===========================================
# Code mode
#===========================================
&set_pp_keywords("Current_Code_Mode" ,"1");
#       Choices
&set_pp_keywords("Normal_Operation_Mode"    ,"0");  #!
&set_pp_keywords("Make_Glass_Mode"          ,"1");  #!
&set_pp_keywords("Shock_Tube_Test_Mode"     ,"2");  #!
&set_pp_keywords("Blastwave_Test_Mode"      ,"3");  #! Sedov-Taylor blastwave test
&set_pp_keywords("Surface_Tension_Test_Mode","4");  #!
&set_pp_keywords("KHI_Test_Mode"            ,"5");  #! Kelvin-Helmholtz instability
&set_pp_keywords("RTI_Test_Mode"            ,"6");  #! Rayleigh-Taylor instability
&set_pp_keywords("Blob_Test_Mode"           ,"7");  #! future
&set_pp_keywords("Gresho_Vortex_Test_Mode"  ,"8");  #! future
&set_pp_keywords("Noh_Problem_Mode"         ,"9");  #! future; see Springel (2010)
&set_pp_keywords("Keplerian_Ring_Test_Mode" ,"10"); #! future
&set_pp_keywords("Evrard_Test_Mode"         ,"11"); #! 
# &set_pp_keywords("Iliev2006_Test1_Mode"       ,"1");  #!
# &set_pp_keywords("Iliev2006_Test2_Mode"       ,"2");  #!
# &set_pp_keywords("Hubber2013_Test1_Mode"      ,"17"); #!
# &set_pp_keywords("Hubber2013_Test2_Mode"      ,"18"); #!
# &set_pp_keywords("Hubber2013_Test3_Mode"      ,"19"); #!
# &set_pp_keywords("Jubelgas2004_Test1_Mode"    ,"20"); #!
# &set_pp_keywords("Jubelgas2004_Test2_Mode"    ,"21"); #!
# [Notes]
#  - Iliev's tests are standard radiation 
#    hydrodynamics tests
#  - Hubber's tests are sink particle tests.
#  - Jubelgas's tests are thermal conduction tests.

#===========================================
#  Base settings
#===========================================
&set_pp_keywords("Restart_Mode"   ,"$restart"); #!
&set_pp_keywords("Make_Rate_Tables" ,"$mktable"); #!

#===========================================
#  Basic setting for simulation 
#===========================================
# Simulation dimension
&set_pp_keywords("Simulation_Dimension","3");
&error_condition("(Simulation_Dimension != 1) && " .
                 "(Simulation_Dimension != 2) && " .
                 "(Simulation_Dimension != 3)");

#==========================
# Hydrodynamics
#==========================
&set_pp_keywords("Coupled_with_Hydro","1");

#==========================
# ISM Model
#==========================
# [1] Dust
&set_pp_keywords("Consider_Size_Distribution","0");
# [2] Property of ISM
&set_pp_keywords("Multicomponent_ISM","1");
# Note that if Multicomponent_ISM=0, gam_fixed and
# mu_fixed are used.
# [3] Gas Temperature
&set_pp_keywords("Bound_Gas_Temperature" ,"1");

&set_pp_keywords("Use_Temperature_Floor" ,"1");
&set_pp_keywords("Temperature_Floor_Type","2");
&set_pp_keywords("Floor_By_EpsG"         ,"1"); #!
&set_pp_keywords("Floor_By_MSPH"         ,"2"); #!

# [4] Treatment of Energy Equation
&set_pp_keywords("Treatment_EnergyEq" ,"2");
&set_pp_keywords("Adiabatic"          ,"1"); #!
&set_pp_keywords("Isothermal"         ,"2"); #!
&set_pp_keywords("Barotropic"         ,"3"); #!
&set_pp_keywords("Cooling"            ,"4"); #!
&set_pp_keywords("Radtr_Cooling"      ,"5"); #!

&set_pp_keywords("Radiation_Force"    ,"1");
# Note that RadiationForce=1 only makes sense if
# Treatment_EnergyEq=Radtr_Cooling.

&set_pp_keywords("Update_Chemical_Abundance","1");
&set_pp_keywords("Update_Gas_Temperature"   ,"1");
&set_pp_keywords("Update_Dust_Temperature"  ,"1");
# Note that these macros only make sense if
# Treatment_EnergyEq=Cooling/Radtr_Cooling.

# [5] Other physical processes
&set_pp_keywords("Thermal_Conduction","0");
&set_pp_keywords("Thermal_Conductivity_Model","1");
# Choice
&set_pp_keywords("Constant_Conductivity","1"); #!
&set_pp_keywords("Spitzer_Conductivity" ,"2"); #!


&error_condition("((Treatment_EnergyEq == Cooling) || " .
                  "(Treatment_EnergyEq == Radtr_Cooling)) && " .
                 "((Update_Chemical_Abundance == 0) && " .
                  "(Update_Gas_Temperature == 0) && " .
                  "(Update_Dust_Temperature == 0))" );
# In this case, calc_chemical_reactions() does not anything!

&error_condition("Update_Chemical_Abundance == 0");
# The current version of calc_chemical_reactions() does not
# support Update_Chemical_Abundance=0, because subcycle timestep
# is determined by chemical timescales only. Therefore, if this
# option = 1, subcycle timestep can become larger than internal
# energy timescale.


&error_condition("(Multicomponent_ISM == 0) && " .
                 "((Treatment_EnergyEq == Cooling) || " .
                 " (Treatment_EnergyEq == Radtr_Cooling))");
# This code does not support radiative cooling for one-component ISM.

#==========================
# Timestepping
#==========================
&set_pp_keywords("Global_Timestepping","1");
&set_pp_keywords("Timescale_Dyn" ,"1"); #!
&set_pp_keywords("Timescale_Engy","2"); #!
&set_pp_keywords("Timescale_Chem","3"); #!
# [Note]
#  Global_Timestepping specifies the timing of 
#  update of hydrodynamic quantities such as velocity.

#==========================
#  Gravity
#==========================
&set_pp_keywords("External_Gravity","0");
&set_pp_keywords("Self_Gravity"    ,"1");

&set_pp_keywords("Force_Error_Analysis","0");

#===================
# SPH formalizm
#===================
# [1] SPH formalizm
&set_pp_keywords("SPH_Formulation"    ,"3");
&set_pp_keywords("Standard_SPH"       ,"1"); #! 
&set_pp_keywords("Springel2005"       ,"2"); #! 
&set_pp_keywords("Hopkins2013_Energy" ,"3"); #!
&set_pp_keywords("Hopkins2013_Entropy","4"); #! future
&set_pp_keywords("Saitoh_Makino_2012" ,"5"); #! future
# [Notes]
#   The meaning of each choice is as follows:
#   Standard_SPH        => classical SPH (not conservative)
#   Springel2005        => Gadget-2; conservative SPH
#   Hopkins2013_Energy  => Pressure-Energy formulation
#   Hopkins2013_Entropy >= Pressure-Entropy formulation
#   Saitoh_Makion_2013  => DISPH

# [1] SPH kernel
&set_pp_keywords("SPH_Kernel"       ,"2");
&set_pp_keywords("M4_Cubic_Spline"  ,"1"); #!
&set_pp_keywords("M4_CS_w_TC92"     ,"2"); #!
&set_pp_keywords("M5_Quartic_Spline","3"); #!
&set_pp_keywords("M6_Quintic_Spline","4"); #!
&set_pp_keywords("Core_Triangle"    ,"5"); #!
&set_pp_keywords("HOCT4"            ,"6"); #! future
&set_pp_keywords("Wendland_C2"      ,"7"); #! future
&set_pp_keywords("Wendland_C4"      ,"8"); #! future
&set_pp_keywords("Wendland_C6"      ,"9"); #! future

# [2] Treatment of artificial viscosity
&set_pp_keywords("Artificial_Viscosity" ,"2");
&set_pp_keywords("Monaghan_Gingold_1983","1"); #!
&set_pp_keywords("Monaghan_1997"        ,"2"); #!
# [3] Swtiches for improving accuracy
&set_pp_keywords("Balsara_Switch"        ,"0");
&set_pp_keywords("Morris_Monaghan_Switch","0");
&set_pp_keywords("Cullen_Dehnen_Switch"  ,"0");
# Error check
&error_condition("(Balsara_Switch == 1) && " .
                 "(Cullen_Dehnen_Switch == 1)");
&error_condition("(Morris_Monaghan_Switch == 1) && " .
                 "(Cullen_Dehnen_Switch == 1)");


#######################################################################
#===============================================
# Configuration for Make_Glass_Mode
#===============================================
# Simulation settings.
#  - Simulation_Type is enforced to be Local_Simulation.
&conditionally_overwrite("Simulation_Dimension","3",
                         "Current_Code_Mode == Make_Glass_Mode");
# Hydrodynamics
#  - Hydrodynamics is solved.
#  - Isothermal is assumed.
&conditionally_overwrite("Coupled_with_Hydro","1",
                         "Current_Code_Mode == Make_Glass_Mode");
&conditionally_overwrite("Treatment_EnergyEq","2",
                         "Current_Code_Mode == Make_Glass_Mode");
# Gravity
#  - No gravity.
&conditionally_overwrite("External_Gravity","0",
                         "Current_Code_Mode == Make_Glass_Mode");
&conditionally_overwrite("Self_Gravity","0",
                         "Current_Code_Mode == Make_Glass_Mode");

#========================================
# Configuration for Shock tube problem
#========================================
# Simulation Type
&conditionally_overwrite("Simulation_Dimension","3",
                         "Current_Code_Mode == Shock_Tube_Test_Mode");
# Hydrodynamics
&conditionally_overwrite("Multicomponent_ISM","0",
                         "Current_Code_Mode == Shock_Tube_Test_Mode");
&conditionally_overwrite("Bound_Gas_Temperature","0",
                         "Current_Code_Mode == Shock_Tube_Test_Mode");
&conditionally_overwrite("Treatment_EnergyEq","1",
                         "Current_Code_Mode == Shock_Tube_Test_Mode");
# Timestepping
&conditionally_overwrite("Global_Timestepping","1",
                         "Current_Code_Mode == Shock_Tube_Test_Mode");
# Gravity
&conditionally_overwrite("External_Gravity","0",
                         "Current_Code_Mode == Shock_Tube_Test_Mode");
&conditionally_overwrite("Self_Gravity","0",
                         "Current_Code_Mode == Shock_Tube_Test_Mode");

# #===================================
# # Configuration for Iliev06 Test1
# #===================================
# # Basic
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# 
# # Hydrodynamics
# &conditionally_overwrite("Coupled_with_Hydro","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","5",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Radiation_Force","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# # ISM
# &conditionally_overwrite("Update_Chemical_Abundance","1",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Update_Gas_Temperature","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Update_Dust_Temperature","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Use_Temperature_Floor","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# 
# # Gravity
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Self_Gravity","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# 
# # Chemical Reactions
# &conditionally_overwrite("Chemical_Reaction1","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction2","1",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction3","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction4","1",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction5","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction6","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction7","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction8","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction9","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Chemical_Reaction10","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# 
# # Radiative Processes
# &conditionally_overwrite("H2_Line_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("H2_CD_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("HI_Bremsstrahlung_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("HI_CI_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("HI_Recombination_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("HI_CE_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("Dust_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("H2_Formation_on_Dust_Heating","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# &conditionally_overwrite("H2_Formation_in_Gas_Heating","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == Iliev2006_Test1_Mode");
# 
# #===================================
# # Configuration for Iliev06 Test2
# #===================================
# # Basic
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# 
# # Radiation Hydrodynamics
# &conditionally_overwrite("Coupled_with_Hydro","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","5",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Radiation_Force","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# 
# # ISM
# &conditionally_overwrite("Update_Chemical_Abundance","1",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Update_Gas_Temperature","1",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Update_Dust_Temperature","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Use_Temperature_Floor","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# 
# # Gravity
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Self_Gravity","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# 
# # Chemical Reactions
# &conditionally_overwrite("Chemical_Reaction1","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction2","1",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction3","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction4","1",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction5","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction6","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction7","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction8","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction9","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Chemical_Reaction10","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# 
# # Radiative Processes
# &conditionally_overwrite("H2_Line_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("H2_CD_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("HI_Bremsstrahlung_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("HI_CI_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("HI_Recombination_Cooling","1",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("HI_CE_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("Dust_Cooling","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("H2_Formation_on_Dust_Heating","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# &conditionally_overwrite("H2_Formation_in_Gas_Heating","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == Iliev2006_Test2_Mode");
# 
# 
# #========================================
# # Configuration for Point explosion Test
# #========================================
# # Basic
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# 
# # Simulation Type
# &conditionally_overwrite("Simulation_Dimension","3",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Simulation_Type","2",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Virtual_Box_Type","1",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("X_Boundary","1",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Y_Boundary","1",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Z_Boundary","1",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# 
# # Hydrodynamics
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","1",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# 
# # Timestepping
# &conditionally_overwrite("Global_Timestepping","1",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# 
# # Gravity
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# &conditionally_overwrite("Self_Gravity","0",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# 
# # Others
# &conditionally_overwrite("Escape_Particle_Monitor","0",
#                          "Current_Code_Mode == Point_Explosion_Test_Mode");
# 
# #========================================
# # Configuration for Surface Tension Test 
# #========================================
# # Basic
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# 
# # Simulation Type
# &conditionally_overwrite("Simulation_Dimension","2",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Simulation_Type","2",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Virtual_Box_Type","1",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("X_Boundary","1",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Y_Boundary","1",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Z_Boundary","2",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# 
# # Hydrodynamics
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","1",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# 
# # Timestepping
# &conditionally_overwrite("Global_Timestepping","1",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# 
# # Gravity
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Self_Gravity","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# &conditionally_overwrite("Escape_Particle_Monitor","0",
#                          "Current_Code_Mode == Surface_Tension_Test_Mode");
# 
# #=======================================================
# # Configuration for Kelvin-Helmholtz Instability Test 
# #=======================================================
# # Basic
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# 
# # Simulation Type
# &conditionally_overwrite("Simulation_Dimension","2",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Simulation_Type","2",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Virtual_Box_Type","1",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("X_Boundary","1",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Y_Boundary","1",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Z_Boundary","2",
#                          "Current_Code_Mode == KHI_Test_Mode");
# 
# # Hydrodynamics
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","1",
#                          "Current_Code_Mode == KHI_Test_Mode");
# 
# # Timestepping
# &conditionally_overwrite("Global_Timestepping","1",
#                          "Current_Code_Mode == KHI_Test_Mode");
# 
# # Gravity
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Self_Gravity","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# &conditionally_overwrite("Escape_Particle_Monitor","0",
#                          "Current_Code_Mode == KHI_Test_Mode");
# 
# #================================================
# # Configuration for 2D-Keplerian Ring Test
# #================================================
# # Basic
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# 
# # Simulation Type
# &conditionally_overwrite("Simulation_Dimension","2",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# &conditionally_overwrite("Simulation_Type","1",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# 
# # Hydrodynamics
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","2",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# 
# # Timestepping
# &conditionally_overwrite("Global_Timestepping","1",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# 
# # Gravity
# &conditionally_overwrite("External_Gravity","1",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# &conditionally_overwrite("Self_Gravity","0",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# &conditionally_overwrite("Escape_Particle_Monitor","0",
#                          "Current_Code_Mode == Keplerian_Ring_Test_Mode");
# 
# #================================================
# # Configuration for Evrard Test Mode
# #================================================
# # Basic
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# 
# # Simulation Type
# &conditionally_overwrite("Simulation_Dimension","3",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# &conditionally_overwrite("Simulation_Type","1",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# 
# # Hydrodynamics
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","1",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# 
# # Timestepping
# &conditionally_overwrite("Global_Timestepping","1",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# 
# # Gravity
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# &conditionally_overwrite("Self_Gravity","1",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# &conditionally_overwrite("Escape_Particle_Monitor","0",
#                          "Current_Code_Mode == Evrard_Test_Mode");
# 
# #================================================
# # Configuration for Isothermal Evrard Test Mode
# #================================================
# # Basic
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# 
# # Simulation Type
# &conditionally_overwrite("Simulation_Dimension","3",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# &conditionally_overwrite("Simulation_Type","1",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# 
# # Hydrodynamics
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","2",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# 
# # Timestepping
# &conditionally_overwrite("Global_Timestepping","1",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# 
# # Gravity
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# &conditionally_overwrite("Self_Gravity","1",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# &conditionally_overwrite("Escape_Particle_Monitor","0",
#                          "Current_Code_Mode == Isothermal_Evrard_Test_Mode");
# 
# 
# 
# #===============================================
# # Configuration for Hubber2013_Test1_Mode
# #===============================================
# # Simulation settins.
# #  - DM,star are OFF.
# #  - sink is ON.
# #  - Global_Simulation is enforced.
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# &conditionally_overwrite("Sink_Component","1",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# &conditionally_overwrite("Simulation_Type","1",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# 
# # Hydrodynamics
# #  - Hydrodynamics solved.
# #  - Multicomponent_ISM ON.
# #  - Isothermal assumed.
# &conditionally_overwrite("Coupled_with_Hydro","1",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# &conditionally_overwrite("Multicomponent_ISM","1",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","2",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# 
# # Gravity
# #  - External gravity OFF.
# #  - Self-gravity ON.
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# &conditionally_overwrite("Self_Gravity","1",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# 
# # Sink properties
# #  - sink creation forbidden.
# #  - sink movement and accretion to sink allowed.
# &conditionally_overwrite("Allow_Sink_Creation","0",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# &conditionally_overwrite("Allow_Sink_Motion","1",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# &conditionally_overwrite("Allow_Sink_Accretion","1",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# 
# # Others
# #  - erase_fluctuation() not called.
# &conditionally_overwrite("Erase_Fluctuaion","0",
#                          "Current_Code_Mode == Hubber2013_Test1_Mode");
# 
# #===============================================
# # Configuration for Hubber2013_Test2_Mode
# #===============================================
# # Simulation settins.
# #  - DM,star are OFF.
# #  - sink is ON.
# #  - Global_Simulation is enforced.
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# &conditionally_overwrite("Sink_Component","1",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# &conditionally_overwrite("Simulation_Type","1",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# 
# # Hydrodynamics
# #  - Hydrodynamics solved.
# #  - Multicomponent_ISM ON.
# #  - Isothermal assumed.
# &conditionally_overwrite("Coupled_with_Hydro","1",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# &conditionally_overwrite("Multicomponent_ISM","1",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","2",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# 
# # Gravity
# #  - External gravity OFF.
# #  - Self-gravity ON.
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# &conditionally_overwrite("Self_Gravity","1",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# 
# # Sink properties
# #  - sink creation allowed.
# #  - sink movement and accretion to sink allowed.
# &conditionally_overwrite("Allow_Sink_Creation","1",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# &conditionally_overwrite("Allow_Sink_Motion","1",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# &conditionally_overwrite("Allow_Sink_Accretion","1",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# 
# # Others
# #  - erase_fluctuation() not called.
# &conditionally_overwrite("Erase_Fluctuaion","0",
#                          "Current_Code_Mode == Hubber2013_Test2_Mode");
# 
# #===============================================
# # Configuration for Hubber2013_Test3_Mode
# #===============================================
# # Simulation settins.
# #  - DM,star are OFF.
# #  - sink is ON.
# #  - Global_Simulation is enforced.
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("Sink_Component","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("Simulation_Type","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# 
# # Hydrodynamics
# #  - Hydrodynamics solved.
# #  - Multicomponent_ISM OFF.
# #  - Isothermal assumed.
# &conditionally_overwrite("Coupled_with_Hydro","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","3",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# 
# # Gravity
# #  - External gravity OFF.
# #  - Self-gravity ON.
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("Self_Gravity","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# 
# # Sink properties
# #  - sink creation allowed.
# #  - sink movement and accretion to sink allowed.
# &conditionally_overwrite("Allow_Sink_Creation","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("Allow_Sink_Motion","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("Allow_Sink_Accretion","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("SSC_Density_Threshold","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("SSC_Jeans_Condition","0",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("SSC_Minimum_Phi","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("SSC_No_Overlap","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# &conditionally_overwrite("SSC_HiLL_Condition","1",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# 
# # Others
# #  - erase_fluctuation() not called.
# &conditionally_overwrite("Erase_Fluctuaion","0",
#                          "Current_Code_Mode == Hubber2013_Test3_Mode");
# 
# #===============================================
# # Configuration for Jubelgas2004_Test1_Mode
# #===============================================
# # Simulation settings.
# #  - DM,star,sink are OFF.
# #  - Simulation_Dimension=3.
# #  - Local_Simulation is enforced.
# #  - Periodic boundary conditions are applied for y,z-directions,
# #    while vacuum boundary is applied for x-direction.
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Simulation_Dimension","3",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Simulation_Type","2",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Virtual_Box_Type","1",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("X_Boundary","2",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Y_Boundary","1",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Z_Boundary","1",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# 
# # Hydrodynamics
# #  - Hydrodynamics not solved.
# #  - Multicomponent_ISM OFF.
# #  - Bound_Gas_Temperature OFF.
# #  - Adiabatic assumed.
# #  - Thermal Conduction ON.
# &conditionally_overwrite("Coupled_with_Hydro","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","1",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Thermal_Conduction","1",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# 
# # Gravity
# #  - External gravity OFF.
# #  - Self-gravity OFF.
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Self_Gravity","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# &conditionally_overwrite("Escape_Particle_Monitor","0",
#                          "Current_Code_Mode == Jubelgas2004_Test1_Mode");
# 
# #===============================================
# # Configuration for Jubelgas2004_Test2_Mode
# #===============================================
# # Simulation settings.
# #  - DM,star,sink are OFF.
# #  - Simulation_Dimension=3.
# #  - Global_Simulation is enforced.
# &conditionally_overwrite("DarkMatter_Component","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Stellar_Component","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Sink_Component","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Simulation_Dimension","3",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Simulation_Type","1",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Virtual_Box_Type","1",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# 
# # Hydrodynamics
# #  - Hydrodynamics not solved.
# #  - Multicomponent_ISM OFF.
# #  - Bound_Gas_Temperature OFF.
# #  - Adiabatic assumed.
# #  - Thermal Conduction ON.
# &conditionally_overwrite("Coupled_with_Hydro","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Multicomponent_ISM","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Bound_Gas_Temperature","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Treatment_EnergyEq","1",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Thermal_Conduction","1",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# 
# # Gravity
# #  - External gravity OFF.
# #  - Self-gravity OFF.
# &conditionally_overwrite("External_Gravity","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Self_Gravity","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# 
# # Others
# &conditionally_overwrite("Erase_Fluctuation","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# &conditionally_overwrite("Escape_Particle_Monitor","0",
#                          "Current_Code_Mode == Jubelgas2004_Test2_Mode");
# 
# #===============================================
# #   Final check for the keywords
# #===============================================
# # [Notes]
# #   - In 1D calculation, we use only x coordinates.
# #     Therefore, y,z boudanry conditions should be vaccum B.C..
# #   - In 2D calculation, we use x,y coordinates only.
# #     Threfore, z B.C. force to be vaccum B.C..
# &conditionally_overwrite("Z_Boundary","2",
#                          "(Simulation_Dimension == 1) || " .
# 			 "(Simulation_Dimension == 2)");
# &conditionally_overwrite("Y_Boundary","2",
#                          "(Simulation_Dimension == 1)");

#-----------------------------------------------------------
# [3] Open file & Output
#-----------------------------------------------------------
open File, ">preprocess_keywords.h" or die "cannot make the file.";
   # [3-1] write simple macro definitions
   if ($pointer ne "-1") {
      for (my $i=0;$i<=$pointer;$i++) {
         my $output = "#define " . $ordering[$i] . " (" . $pp_keywords{$ordering[$i]} . ")";
         print File "$output\n";
      }
   }
   print File "\n\n";
   # [3-2] write #pragma once
   print File "#pragma once\n";
   # [3-3] write conditionally-defined macros
   print File "//=============================\n";
   print File "// Conditionally-defined macros\n";
   print File "//=============================\n";
   if ($pointer_wc ne "-1") {
      for (my $i=0;$i<=$pointer_wc;$i++) {
         my $output = "#if (" . $pp_conditions_wc[$i] . ")\n"
	 	    . "#ifdef " . $pp_keywords_wc[$i] . "\n"
		    . "#undef " . $pp_keywords_wc[$i] . "\n"
		    . "#endif\n"
                    . "#define " . $pp_keywords_wc[$i] . " (" . $pp_values_wc[$i] . ")\n"
                    . "#endif";
         print File "$output\n\n";
      }
   }
   print File "\n\n";
   # [3-4] write error detect script
   print File "//=============================\n";
   print File "// Error conditions            \n";
   print File "//=============================\n";
   if ($pointer_err ne "-1") {
      for (my $i=0;$i<=$pointer_err;$i++) {
         my $output = "#if (" . $err_list[$i] . ")\n"
                    . "#error\n"
                    . "#endif";
         print File "$output\n\n";
      }
   }
close File;


sub set_pp_keywords (\$\$) {
  my $key   = shift;
  my $value = shift;

  $pointer++;
  $pp_keywords{$key} = $value;
  $ordering[$pointer] = $key;

}

# Conditionally overwrite
sub conditionally_overwrite (\$\$\$) {
  my $key       = shift;
  my $value     = shift;
  my $condition = shift;

  $pointer_wc++;
  $pp_keywords_wc[$pointer_wc]   = $key;
  $pp_values_wc[$pointer_wc]     = $value;
  $pp_conditions_wc[$pointer_wc] = $condition;

}

# Error Check
sub error_condition (\$) {
  my $err_condition = shift;

  $pointer_err++;
  $err_list[$pointer_err] = $err_condition;

}
