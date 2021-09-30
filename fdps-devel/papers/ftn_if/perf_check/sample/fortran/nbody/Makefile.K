#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
FDPS_LOC = ../../../../../../../../
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_FTN_MOD_DIR = $(FDPS_LOC)/src/fortran_interface/modules
FDPS_FTN_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_ftn_if.py

# (ii) Variables to specify compilers and compile options
# Serial or OpenMP cases
FC=time frtpx
CXX=time FCCpx
# MPI case
#FC=time -p mpifrtpx
#CXX=time -p mpiFCCpx
# [Option 1] w/ debug
#FCFLAGS = -X 03 -O0 -g -Haefosux -Free -fs
#CXXFLAGS = -O0 -g -NRtrap -Nquickdbg=subchk -Nquickdbg=heapchk -Nquickdbg=inf_detail $(FDPS_INC) 
# [Option 2] w/ optimization 
FCFLAGS = -X 03 -Kfast -Free -fs
#CXXFLAGS = -g -Kfast,nons,noeval,nofp_contract,nofp_relaxed $(FDPS_INC)
CXXFLAGS = -O3 $(FDPS_INC)
LDFLAGS = --linkfortran
# OpenMP options
#FCFLAGS  += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -openmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -openmp
# MPI options
#FCFLAGS  += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
# Other options
FCFLAGS += -DINTERACTION_FUNCTION_IS_EMPTY


# (iii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K-computer.sh
#WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_10
#WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_11
#WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_12
#WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_13
#WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_14
#WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_15
#WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_16
#WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_17
WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran/work_2_18


#-----------------------------------------------------------------------
#   Source files
#-----------------------------------------------------------------------
%.o : %.F90
	$(FC) $(FCFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

SRC_USER_DEFINED_TYPE = user_defined.F90
SRC_USER = f_main.F90
SRC_FDPS_MOD = $(wildcard $(FDPS_FTN_MOD_DIR)/*.F90)
SRC_FTN = $(SRC_FDPS_MOD) \
	  $(SRC_USER_DEFINED_TYPE) \
	  FDPS_module.F90 \
	  $(SRC_USER)
SRC_CXX = FDPS_ftn_if.cpp \
	  FDPS_Manipulators.cpp \
	  main.cpp

OBJ_USER_DEFINED_TYPE	= $(SRC_USER_DEFINED_TYPE:F90=o)
OBJ_USER		= $(SRC_USER:F90=o)
OBJ_FDPS_MOD		= $(notdir $(SRC_FDPS_MOD:F90=o))
OBJ_FTN			= $(notdir $(SRC_FTN:F90=o))
OBJ_CXX			= $(SRC_CXX:cpp=o)
OBJ			= $(OBJ_FTN) $(OBJ_CXX)

VPATH = $(FDPS_FTN_MOD_DIR)
TARGET = nbody.out

$(TARGET): $(OBJ) result
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)
	mkdir -p $(WORKDIR)
	mv $(TARGET) $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME)  

result:
	mkdir result

$(SRC_CXX) FDPS_module.F90: $(SRC_USER_DEFINED_TYPE)
	$(FDPS_FTN_IF_GENERATOR) $(SRC_USER_DEFINED_TYPE) --output ./ 

FDPS_super_particle.o: FDPS_vector.o FDPS_matrix.o

$(OBJ_USER_DEFINED_TYPE): $(OBJ_FDPS_MOD)

FDPS_module.o: $(OBJ_USER_DEFINED_TYPE)

$(OBJ_USER): $(OBJ_USER_DEFINED_TYPE) FDPS_module.o

clean:
	rm -f *.o *.s *.mod $(TARGET) *.dat

distclean: clean
	rm -f $(SRC_CXX) FDPS_Manipulators.h  FDPS_module.F90 user_defined.hpp
	rm -f *__genmod.f90
	rm -rf result
