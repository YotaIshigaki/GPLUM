#-----------------------------------------------------------------------
#   Compiler and complile options                                       
#-----------------------------------------------------------------------
#HOST = mac
#HOST = armagnac0
HOST = knl

# [1] Mac Env.
ifeq ($(HOST),mac)
FC=gfortran-mp-5
CC=gcc
CXX=g++-mp-5

FDPS_LOC = /Users/name/fdps
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_LIB = 

FCFLAGS = -O0 -Wall
CFLAGS =
CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
LDFLAGS = -lgfortran $(FDPS_LIB)
endif

# [2] armagnac0 & jenever0
ifeq ($(HOST),armagnac0)
# Serial or OpenMP only
FC=gfortran
CC=gcc
CXX=g++
# MPI 
#FC=mpif90
#CC=mpicc
#CXX=mpic++

MPI_LOC  = /usr/lib64/openmpi
MPI_INC  = -I$(MPI_LOC)/lib
MPI_LIB  = -L$(MPI_LOC)/lib -lmpi -lmpi_f90

FDPS_LOC = /nfshome/namekata/codes/fdps
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_LIB = 

# [Option 1] w/o optimization
#FCFLAGS = -std=f2003 -O0 -Wall
#CFLAGS =
#CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
FCFLAGS = -std=f2003 -O3 -ffast-math -funroll-loops -finline-functions
CFLAGS =
CXXFLAGS = -O3 -ffast-math -funroll-loops $(FDPS_INC)
# --- Common part ---
#FCFLAGS   += -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
LDFLAGS = -lgfortran $(MPI_LIB) $(FDPS_LIB)
endif

# [3] knl
ifeq ($(HOST),knl)
# Serial or OpenMP only
#FC=/opt/intel/bin/ifort
#CC=/opt/intel/bin/icc
#CXX=/opt/intel/bin/icpc
# MPI
FC=time -p mpiifort
CC=time -p mpiicc
CXX=time -p mpiicpc

MPI_LOC  = /opt/intel/compilers_and_libraries/linux/mpi
MPI_BIN  = $(MPI_LOC)/bin64
MPI_INC  = -I$(MPI_LOC)/include64
MPI_LIB  = -L$(MPI_LOC)/lib64

FDPS_LOC = /home/namekata/codes/fdps
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_LIB = 

# [Option 1] w/o optimization
#FCFLAGS = -fast -I/opt/intel/include
#CFLAGS = -fast
#CXXFLAGS = -fast -I/opt/intel/include $(FDPS_INC) 
# [Option 2] w/ debug
#FCFLAGS = -O0 -check all -fpe0 -gen_interfaces -ftrapuv -traceback -g
#CFLAGS =
#CXXFLAGS = -ftrapuv -debug all -g -I/opt/intel/include $(FDPS_INC) 
# [Option 3] w/ optimization 
#FCFLAGS = -fast -ipo 
#CFLAGS = -fast -ipo
#CXXFLAGS = -fast -ipo $(FDPS_INC)
FCFLAGS = -O3 -ipo -xMIC-AVX512
CFLAGS = -O3 -ipo -xMIC-AVX512
CXXFLAGS = -O3 -ipo -xMIC-AVX512 $(FDPS_INC)
# --- Common part ---
FCFLAGS  += -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DPARTICLE_SIMULATOR_VECTOR_RANGE_CHECK
LDFLAGS = -L/opt/intel/lib/intel64 -L/usr/lib64 -lifport -lifcore -limf -lsvml -lm -lipgo -lirc -lirc_s $(MPI_LIB) $(FDPS_LIB)
endif

#-----------------------------------------------------------------------
#   Source files                                                        
#-----------------------------------------------------------------------
SRC_FTN = FDPS_vector.F90 \
	  FDPS_matrix.F90 \
	  FDPS_super_particle.F90 \
	  FDPS_time_profile.F90 \
	  user_defined.F90 \
	  FDPS_module.F90 \
	  f_main.F90
SRC_C   =
SRC_CXX = FDPS_ftn_ifc.cpp \
	  FDPS_Manipulators.cpp \
	  main.cpp
SRC     = $(SRC_FTN) $(SRC_C) $(SRC_CXX)

OBJ_FTN = $(addsuffix .o, $(basename $(SRC_FTN))) # $(SRC1:.F90=.o)
OBJ_C   = $(addsuffix .o, $(basename $(SRC_C)))   # $(SRC1:.c=.o)
OBJ_CXX = $(addsuffix .o, $(basename $(SRC_CXX))) # $(SRC3:.cpp=.o)
OBJ     = $(OBJ_FTN) $(OBJ_C) $(OBJ_CXX)

TARGET = a.out

$(TARGET): Makefile $(SRC)
	make obj

obj: Makefile $(SRC)
	make obj_ftn
	#make obj_c
	make obj_cxx
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)
	strip $(TARGET)


obj_ftn: $(SRC_FTN) Makefile
	$(FC) $(FCFLAGS) -c $(SRC_FTN)

obj_c: $(SRC_C) Makefile
	$(CC) $(CFLAGS) -c $(SRC_C)

obj_cxx: $(SRC_CXX) Makefile
	$(CXX) $(CXXFLAGS) -c $(SRC_CXX)

clean:
	rm -f *.o *.s *.mod $(TARGET) *.dat

