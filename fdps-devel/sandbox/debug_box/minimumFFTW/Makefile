#--------------------------------------------------------------------
#  Makefile for "aTSEVB_FDPS.x"
#    2015/12/24  K.Kawai Tokumasu-Lab.
#--------------------------------------------------------------------

#--- PATH for FDPS src directry
PS_DIR   = ../../..
PS_PATH  = -I$(PS_DIR)/src/
PS_PATH += -I$(PS_DIR)/src/particle_mesh/
LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm.a
#LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm_debug.a

#--- add FFTw library
#FFTW_DIR      = /opt/local
FFTW_DIR      = /home/masaki/opt
INCLUDE_FFTW  = -I$(FFTW_DIR)/include/
LIB_FFTW      = -L$(FFTW_DIR)/lib/ -lfftw3f_mpi
LIB_FFTW_STATIC  =   $(FFTW_DIR)/lib/libfftw3f_mpi.a
LIB_FFTW_STATIC +=   $(FFTW_DIR)/lib/libfftw3f.a

#--- compiler settings
#CC = g++
#CC = mpicxx-openmpi-gcc6
CC = mpicxx
CFLAGS = -O3 -ffast-math -funroll-loops -lm -g

FC = gfortran
FFLAGS = -std=2003

#--- enable parallelization in FDPS
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#--- enable user-added debug function
CFLAGS += -DDINFO_TEST
CFLAGS += -DMASK_TEST
CFLAGS += -DVIRIAL_TEST


#--- Fortran compile target (till the objects)
#SRC_F   = (add *.F90 file here)
OBJ_F   = $(SRC_F90:%.f90=%.o)

#--- C++ compile target
PROGRAM = test.out
SRC_CPP = test.cpp
OBJ_CPP = $(SRC_CPP:%.cpp=%.o)
OBJS    = $(OBJ_CPP) $(OBJ_F)

#--- dependency check
DEPS = $(SRC_CPP:%.cpp=%.d)


all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(PS_PATH) $(LIB_FFTW) $(CFLAGS) -o $@ $< $(LIB_PM) $(LIB_FFTW_STATIC)

$(OBJ_CPP): $(SRC_CPP)
	$(CC) -c $(PS_PATH) $(INCLUDE_FFTW) $(CFLAGS) -MMD $<

$(OBJ_F): $(SRC_F)
	$(FC) -c $(FFLAGS) $<
    
static: $(OBJS)
	$(CC) $(PS_PATH) $(CFLAGS) -o $(PROGRAM) $< $(LIB_PM) $(LIB_FFTW_STATIC)
    
.PHONY: clean
clean:
	rm $(PROGRAM) $(DEPS) $(OBJS)


-include $(DEPS)
    
