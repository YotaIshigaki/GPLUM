#PS_PATH = ../../../../../../src/
#PS_PATH = ../../../../my_trunk/src_old04/
PS_PATH = ../../../../my_trunk/src/

INC = -I$(PS_PATH)

#CC = time g++
#CC = time CC
CC = mpicxx
CFLAGS = -O3
CFLAGS += -Wall
CFLAGS += -ffast-math
CFLAGS += -funroll-loops
CFLAGS += -std=c++11
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -DMONAR
#CFLAGS += -DREUSE_LIST_MODE
CFLAGS += -DSHRINK_SIZE_OF_SPJ_SORTED
#CFLAGS += -DCHECK_ENERGY

use_phantom_grape_x86 = yes
#use_gpu_cuda = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

#all:nbody.out nbody.new.out
all:nbody.out

ifeq ($(use_phantom_grape_x86),yes)
PHANTOM_PATH = ../../../../../../src/
INC += -I$(PHANTOM_PATH)/phantom_grape_x86/G5/newton/libpg5/
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CLIBS   = -L$(PHANTOM_PATH)/phantom_grape_x86/G5/newton/libpg5/ -lpg5
endif

ifeq ($(use_gpu_cuda),yes)
CUDA_HOME = /usr/local/cuda
#CUDA_HOME = /gwfefs/opt/x86_64/cuda/7.5
NVCC = time $(CUDA_HOME)/bin/nvcc -Xcompiler="-O3 -Wall -fopenmp -march=core-avx2 -ffast-math -funroll-loops"
INC  += -I$(CUDA_HOME)/samples/common/inc/
CFLAGS  += -DENABLE_GPU_CUDA


#CFLAGS  += -DMULTI_WALK -DSHRINK_SIZE_OF_SPJ_SORTED
#CFLAGS  += -DMULTI_WALK_INDEX -DSHRINK_SIZE_OF_SPJ_SORTED
#CFLAGS  += -DMULTI_WALK_INDEX -DSHRINK_SIZE_OF_SPJ_SORTED -DREUSE_LIST_MODE

#CFLAGS  += -DMULTI_WALK_STREAM -DSHRINK_SIZE_OF_SPJ_SORTED
CFLAGS  += -DMULTI_WALK_INDEX_STREAM -DSHRINK_SIZE_OF_SPJ_SORTED -DREUSE_LIST_MODE 

CLIBS = -L$(CUDA_HOME)/lib64 -lcudart -lgomp
#CUDAFLAGS += -DSANITY_CHECK_CUDA_POINTER
CUDAFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CUDAFLAGS += -DUSE_CUDA_TIMER
#CUDAFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
force_gpu_cuda.o:force_gpu_cuda.cu
	$(NVCC) $(INC) $(CUDAFLAGS) -c -o $@ $< 
OBJS = force_gpu_cuda.o
endif

nbody.out:nbody.cpp $(OBJS)
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(CLIBS)

#nbody.new.out:nbody.cpp $(OBJS)
#	$(CC) $(INC) $(CFLAGS) -DSHRINK_SIZE_OF_SPJ_SORTED -o $@ $^ $(CLIBS)

clean:
	rm -f nbody.out

distclean:
	rm -f nbody.out *.o *~

test:
	# This command is only for FDPS developers.
	./test.py

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
