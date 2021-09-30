PS_PATH = ../../../../../../../../src/
INC = -I$(PS_PATH)

CC=time FCCpx
#CC=time mpiFCCpx
#CFLAGS = -Kfast,nons,noeval,nofp_contract,nofp_relaxed -Ksimd=2 -Krestp=all
CFLAGS = -Kfast -Ksimd=2 -Krestp=all
#CFLAGS += -O0 -g -NRtrap -Nquickdbg=subchk -Nquickdbg=heapchk -Nquickdbg=inf_detail
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -openmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += -DINTERACTION_FUNCTION_IS_EMPTY

# Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/ftn_if/nbody/c++


#use_phantom_grape_x86 = yes
#use_gpu_cuda = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

all:nbody.out

ifeq ($(use_phantom_grape_x86),yes)
INC += -I$(PS_PATH)/phantom_grape_x86/G5/newton/libpg5/
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CLIBS   = -L$(PS_PATH)/phantom_grape_x86/G5/newton/libpg5/ -lpg5
endif

ifeq ($(use_gpu_cuda),yes)
CUDA_HOME = /usr/local/cuda
#CUDA_HOME = /gwfefs/opt/x86_64/cuda/7.5
NVCC = time $(CUDA_HOME)/bin/nvcc -Xcompiler="-O3"
INC  += -I$(CUDA_HOME)/samples/common/inc/
CFLAGS  += -DENABLE_GPU_CUDA
CLIBS = -L$(CUDA_HOME)/lib64 -lcudart -lgomp
force_gpu_cuda.o:force_gpu_cuda.cu
	$(NVCC) $(INC) -c -o $@ $<
OBJS = force_gpu_cuda.o
endif

nbody.out:nbody.cpp $(OBJS)
	#$(CC) $(INC) $(CFLAGS) -Koptmsg=2 -S -o nbody.s $^ $(CLIBS) > compile_log.txt 2>&1
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(CLIBS)
	mkdir -p $(WORKDIR)
	mv $@ $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME)

clean:
	rm -f nbody.out

distclean:
	rm -f nbody.out *.o *~

test: 
	# This command is only for FDPS developers.
	./test.py

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
