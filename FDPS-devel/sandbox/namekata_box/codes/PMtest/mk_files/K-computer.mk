#-----------------------------------------------------------------------
#   Compiler and complile options
#-----------------------------------------------------------------------
FC=mpifrtpx
CC=mpifccpx
CXX=mpiFCCpx

# Library informations
LIBDIR          = /home/ra000008/a03339/fujitsu
LIBpDIR         = $(dir $(LIBDIR))
LIBCOMP         = $(notdir $(LIBDIR)).tar.gz

GSL_INC		= -I/home/ra000008/a03339/fujitsu/gsl-2.1/include
GSL_LIBS	= -L/home/ra000008/a03339/fujitsu/gsl-2.1/lib -lgsl -lgslcblas
GSL_LD_PATH	= ./fujitsu/gsl-2.1/lib
QD_INC		= -I/home/ra000008/a03339/fujitsu/include
QD_LIBS		= -L/home/ra000008/a03339/fujitsu/lib -lqd -lm
QD_LD_PATH      = ./fujitsu/lib
# ~/fujitsu/bin/qd-config --cxxflags
# ~/fujitsu/bin/qd-config --libs
EIGEN_INC	= -I/home/ra000008/a03339/fujitsu/include/eigen3
EIGEN_LIBS	=
FFTW_INC        = -I/home/apps/fftw/3.3.3_1.2.0-19/include
FFTW_LIBS	= -L/home/apps/fftw/3.3.3_1.2.0-19/lib64 -lfftw3_mpi -lfftw3_omp -lfftw3 -lm
FDPS_INC	= -I../../../../src -I../../../../src/particle_mesh
FDPS_LIBS	= -L../../../../src/particle_mesh -lpm

#DEBUG = -g -Nquickdbg
DEBUG = 
CXXFLAGS = -Xg -Kfast,nons,noeval,nofp_contract,nofp_relaxed $(GSL_INC) $(QD_INC) $(EIGEN_INC) $(FDPS_INC) $(FFTW_INC)
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += $(DEBUG)
LDFLAGS = $(QD_LIBS) $(GSL_LIBS) $(EIGEN_LIBS) $(FDPS_LIBS) $(FFTW_LIBS)

#======================================================================
#   Job scheduling command
#======================================================================
QSUB = pjsub --step

#----------------------------------------------------------------------
#   WORKDIR
#----------------------------------------------------------------------
# Normal queue
WORKDIR = /data/ra000008_a/namekata/Debug
#WORKDIR = /data/ra000008_a/namekata/Debug2
#WORKDIR = /data/ra000008_a/namekata/Debug3

# micro
#WORKDIR = /scratch/ra000008/namekata/Debug

#-----------------------------------------------------------------------
#   JOB_FILE
#-----------------------------------------------------------------------
BATCH_DIR                       = ./batch_files/K-computer
DEBUG_FILE_NAME                 = debug.sh
DEBUG_RESTART_FILE_NAME         = debug_restart.sh
MKTABLE_FILE_NAME               = mktable.sh
JOB_FILE_NAME                   = job.sh
RESTART_FILE_NAME               = restart.sh
DEBUG_FILE              = $(BATCH_DIR)/$(DEBUG_FILE_NAME)
DEBUG_RESTART_FILE      = $(BATCH_DIR)/$(DEBUG_RESTART_FILE_NAME)
MKTABLE_FILE            = $(BATCH_DIR)/$(MKTABLE_FILE_NAME)
JOB_FILE                = $(BATCH_DIR)/$(JOB_FILE_NAME)
RESTART_FILE            = $(BATCH_DIR)/$(RESTART_FILE_NAME)


