#--------------------------------------------------------------------
#  Makefile for "aTSEVB_FDPS.x"
#    2015/12/24  K.Kawai Tokumasu-Lab.
#--------------------------------------------------------------------

#--- PATH for FDPS src directry
#------ for paolo in Tokumasu Lab.
#PS_DIR   = /home/kawai/local/FDPS_2.0_20160302
#PS_DIR   = /home/kawai/local/FDPS_2.0_20160125
#PS_DIR   = /home/kawai/local/FDPS_1.1_20150828
PS_DIR 	= /nfshome/namekata/codes/fdps

PS_PATH  = -I$(PS_DIR)/src/
PS_PATH += -I$(PS_DIR)/src/particle_mesh/
#LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm.a
LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm_debug.a

#--- add FFTw library
FFTW_DIR      = /nfshome/namekata/gnu
INCLUDE_FFTW  = -I$(FFTW_DIR)/include/
LIB_FFTW      = -L$(FFTW_DIR)/lib/ -lfftw3f_mpi -lfftw3f
LIB_FFTW_STATIC  =   $(FFTW_DIR)/lib/libfftw3f_mpi.a
LIB_FFTW_STATIC +=   $(FFTW_DIR)/lib/libfftw3f.a

#--- add MD extention
INCLUDE_MD_EXT = -I./md_ext_include/

#--- compiler settings
#CC = g++
CC = mpicxx
#CFLAGS = -O3 -ffast-math -funroll-loops -lm -std=c++11 -g3
CFLAGS = -O0 -Wall -Wextra -lm -std=c++11 -g3 
#CFLAGS = -O0 -Wall -Wextra -lm -std=c++11 -g3 -DPARTICLE_SIMULATOR_DEBUG_PRINT

#--- enable parallelization in FDPS
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#--- enable user-added debug function
CFLAGS += -DDINFO_TEST
CFLAGS += -DMASK_TEST
#CFLAGS += -DVIRIAL_TEST

#--- C++ compile target
PROGRAM = test.out
SRC_CPP = atsevb_main.cpp
#SRC_CPP = g_p3m_test.cpp
OBJ_CPP = $(SRC_CPP:%.cpp=%.o)
OBJS    = $(OBJ_CPP)

#--- dependency check
DEPS = $(SRC_CPP:%.cpp=%.d)


all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(PS_PATH) $(LIB_FFTW) $(CFLAGS) -o $@ $< $(LIB_PM)

$(OBJ_CPP): $(SRC_CPP)
	$(CC) -c $(PS_PATH) $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(CFLAGS) -MMD $<

static: $(OBJS)
	$(CC) $(PS_PATH) $(CFLAGS) -o $(PROGRAM) $< $(LIB_PM) $(LIB_FFTW_STATIC)
    
.PHONY: clean
clean:
	rm $(PROGRAM) $(DEPS) $(OBJS)


-include $(DEPS)
    
