#======================================================================
#   Numerical Libraries and Compilers
#======================================================================
QD_INC = -I/nfshome/namekata/gnu/include
QD_LIB = -L/nfshome/namekata/gnu/lib -lqd -lm -lrt
# /nfshome/namekata/gnu/bin/qd-config --cxxflags
# /nfshome/namekata/gnu/bin/qd-config --libs
GSL_INC = -I/nfshome/namekata/gnu/gsl-2.1/include
GSL_LIB = -L/nfshome/namekata/gnu/gsl-2.1/lib -lgsl -lgslcblas
EIGEN_INC = -I/nfshome/namekata/gnu/include/eigen3
EIGEN_LIB = 
FDPS_INC = -I../../../../src
FDPS_LIB =
LDFLAGS = $(QD_LIB) $(GSL_LIB) $(EIGEN_LIB) $(FDPS_LIB)

#CFLAGS_COMMON = -std=c++11 -O3 -ffast-math -funroll-loops $(QD_INC) $(GSL_INC) $(EIGEN_INC) $(FDPS_INC)
CFLAGS_COMMON = -std=c++11 -Wall -Wextra -ftrapv -fexceptions -g3 $(QD_INC) $(GSL_INC) $(EIGEN_INC) $(FDPS_INC)
# [1] Serial
#CC = g++
#CFLAGS = $(CFLAGS_COMMON) 
#CFLAGS = $(CFLAGS_COMMON) -DPARTICLE_SIMULATOR_DEBUG_PRINT
# [2] OpenMP
#CC = g++
#CFLAGS = $(CFLAGS_COMMON) -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# [3] MPI
CC = mpicxx
CFLAGS = $(CFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL 
# [4] MPI + OpenMP
#CC = mpicxx
#CFLAGS = $(CFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

#======================================================================
#   Job scheduling command
#======================================================================
QSUB = pjsub --step

#----------------------------------------------------------------------
#   WORKDIR
#----------------------------------------------------------------------
WORKDIR = /data/ra000008/a03339/Debug
#WORKDIR = /data/ra000008_3/a03339/Debug2
#WORKDIR = /data/ra000008_3/a03339/Debug3

#-----------------------------------------------------------------------
#   JOB_FILE
#-----------------------------------------------------------------------
BATCH_DIR                       = ./batch_files/K-computer
DEBUG_FILE_NAME                 = debug.sh
DEBUG_RESTART_FILE_NAME         = debug_restart.sh
MKTABLE_FILE_NAME               = mktable.sh
JOB_FILE_NAME                   = job.sh
RESTART_FILE_NAME               = restart.sh
DEBUG_FILE              = ./$(BATCH_DIR)/$(DEBUG_FILE_NAME)
DEBUG_RESTART_FILE      = ./$(BATCH_DIR)/$(DEBUG_RESTART_FILE_NAME)
MKTABLE_FILE            = ./$(BATCH_DIR)/$(MKTABLE_FILE_NAME)
JOB_FILE                = ./$(BATCH_DIR)/$(JOB_FILE_NAME)
RESTART_FILE            = ./$(BATCH_DIR)/$(RESTART_FILE_NAME)

