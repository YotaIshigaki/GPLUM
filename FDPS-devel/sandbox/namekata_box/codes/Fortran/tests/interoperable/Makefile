#-----------------------------------------------------------------------
#   Compiler and complile options                                       
#-----------------------------------------------------------------------
HOST = armagnac0
#HOST = knl

# armagnac0 & jenever0
ifeq ($(HOST),armagnac0)
# Serial or OpenMP only
FC=gfortran
CXX=g++

FDPS_LOC = ~/fdps
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_LIB = 
FDPS_FTN_MOD_DIR = ~/fdps/src/fortran_interface/modules

# [Option 1] w/o optimization
#FCFLAGS = -std=f2003 -O0 -Wall
#CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
FCFLAGS = -std=f2003 -O3 -ffast-math -funroll-loops -finline-functions
CXXFLAGS = -O3 -ffast-math -funroll-loops $(FDPS_INC)
# --- Common part ---
LDFLAGS = -lgfortran $(FDPS_LIB)
endif

# knl
ifeq ($(HOST),knl)
# Serial or OpenMP only
FC=/opt/intel/bin/ifort
CXX=/opt/intel/bin/icpc

FDPS_LOC = ~/fdps
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_LIB = 
FDPS_FTN_MOD_DIR = ~/fdps/src/fortran_interface/modules

# [Option 1] w/o optimization
#FCFLAGS = -fast -I/opt/intel/include
#CXXFLAGS = -fast -I/opt/intel/include $(FDPS_INC) 
# [Option 2] w/ debug
#FCFLAGS = -O0 -check all -fpe0 -gen_interfaces -ftrapuv -traceback -g
#CXXFLAGS = -ftrapuv -debug all -g -I/opt/intel/include $(FDPS_INC) 
# [Option 3] w/ optimization 
#FCFLAGS = -fast -ipo 
#CXXFLAGS = -fast -ipo $(FDPS_INC)
FCFLAGS = -O3 -ipo -xMIC-AVX512
CXXFLAGS = -O3 -ipo -xMIC-AVX512 $(FDPS_INC)
# --- Common part ---
LDFLAGS = -L/opt/intel/lib/intel64 -L/usr/lib64 -lifport -lifcore -limf -lsvml -lm -lipgo -lirc -lirc_s $(MPI_LIB) $(FDPS_LIB)
endif

#-----------------------------------------------------------------------
#   Source files                                                        
#-----------------------------------------------------------------------
SRC =   $(FDPS_FTN_MOD_DIR)/FDPS_vector.F90 \
	$(FDPS_FTN_MOD_DIR)/FDPS_matrix.F90 \
	$(FDPS_FTN_MOD_DIR)/FDPS_super_particle.F90 \
	$(FDPS_FTN_MOD_DIR)/FDPS_time_profile.F90 \
	user_defined.F90

OBJ = $(addsuffix .o, $(basename $(notdir $(SRC))))


all: Makefile user_defined.F90
	$(FC) $(FCFLAGS) -c $(SRC)

clean:
	rm -f *.o *.s *.mod $(TARGET) *.dat

