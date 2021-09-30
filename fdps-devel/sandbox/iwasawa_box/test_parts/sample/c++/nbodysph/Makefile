#PS_PATH = -I ../../../src/
PS_PATH = -I../../../../my_trunk/src/
INC = $(PS_PATH)

#CC = time g++
CC = time mpicxx
#CC = mpic++-mpich-gcc49
CFLAGS = -O3
CFLAGS += -ffast-math
CFLAGS += -funroll-loops 
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CFLAGS += -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -DREUSE_LIST_MODE
#CFLAGS += #-Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith

#use_gpu_cuda = yes

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = sph.out

ifeq ($(use_gpu_cuda),yes)
CPPOBJS += force_gpu.o
CUDA_HOME = /usr/local/cuda
NVCC = time $(CUDA_HOME)/bin/nvcc -Xcompiler="-O3"
INC  += -I$(CUDA_HOME)/samples/common/inc/
CFLAGS  += -DENABLE_GPU_CUDA
#CFLAGS  += -DMULTI_WALK
CFLAGS  += -DMULTI_WALK_INDEX
LIBS = -L$(CUDA_HOME)/lib64 -lcudart -lgomp
endif

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH) $(INC)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	$(CC) -c $< $(CFLAGS) $(WARNINGS) $(PS_PATH) $(INC)
	@echo "[$< OK]"

clean:
	-rm *.out *.o


ifeq ($(use_gpu_cuda),yes)
force_gpu.o:force_gpu.cu
	$(NVCC) $(INC) -c -o $@ $<

endif
