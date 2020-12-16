PZCL_KERNEL_DIRS=kernel.sc32

TARGET=nbody_pezy.out

CPPSRC = nbody.cpp

#INC_DIR = -I/home/masaki/project/fdps/src/ 
INC_DIR = -I../../../../src/
INC_DIR += -I/usr/mpi/gcc/openmpi-1.8.4/include/

LIB_DIR =  -lm
LIB_DIR += -L/usr/mpi/gcc/openmpi-1.8.4/lib64

DEFS =  -DENABLE_PEZY 
DEFS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#DEFS += -DMONAR
DEFS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
OMP += -fopenmp

PZCL_STATIC_LIB ?=1
OBJ_DIR ?= obj

#PZSDK_PATH ?= /opt/pzsdk.ver2.1
PZSDK_PATH ?= /opt/pzsdk.ver3.0
PZSDK_INC_PATH ?= $(PZSDK_PATH)/inc
PZSDK_LIB_PATH ?= $(PZSDK_PATH)/lib

INC_DIR += -I$(PZSDK_INC_PATH)

ifeq ($(PZCL_STATIC_LIB),1)
PZCL_LIB ?= $(PZSDK_LIB_PATH)/libpzcl.a
else
PZCL_LIB ?= -lpzcl
endif

LIB_DIR += -L$(PZSDK_LIB_PATH)
LIBS    += $(PZCL_LIB) -lpthread -ldl -lstdc++

PZCL_KERNEL_DIRS ?= 

#CC     ?= g++
CC     = mpicxx

CCOPT  ?= -D__LINUX__ -O2
CCOPT  += -march=core-avx2
CCOPT  += -DFAST_WALK_AVX2
CCOPT  += -Wall
#CCOPT += -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE

###########################################
# Preprocess
###########################################
CPPOBJS= ${patsubst %.cpp, ${OBJ_DIR}/%.o, ${notdir ${CPPSRC}}}
OBJS   = ${CPPOBJS}
DEPS   = ${OBJS:.o=.d} 

all: make_pzcl_kernel ${OBJ_DIR} ${OBJS} ${TARGET} 

nbody.s:nbody.cpp
	${CC} $(OMP) ${DEFS} ${CCOPT} ${INC_DIR} -c -S $< -o $@

clean : clean_pzcl_kernel
	-rm -f ${TARGET} ${OBJS} $(DEPS)
	-rm -rf ${OBJ_DIR}

${TARGET}: ${OBJS}
	${CC} $(OMP) ${DEFS} -o ${TARGET} ${OBJS} ${LIB_DIR} ${LIBS} ${LDOPT}

$(CPPOBJS) : ${OBJ_DIR}/%.o : %.cpp 
	${CC} $(OMP) ${DEFS} ${CCOPT} ${INC_DIR} -c $< -o $@
	${CC} $(OMP) ${DEFS} ${CCOPT} ${INC_DIR} -MM -MT$@ $< > ${OBJ_DIR}/$*.d 

${OBJ_DIR}:
	@if [ ! -d ${OBJ_DIR} ]; then \
		echo "mkdir ${OBJ_DIR}"; mkdir ${OBJ_DIR}; \
	fi

make_pzcl_kernel : $(PZCL_KERNEL_DIRS)
	@for d in $(PZCL_KERNEL_DIRS); \
	do \
	(cd $$d && \
	echo "" && \
	echo "*******************************************************" && \
	echo "* making in ./$$d" && \
	echo "*******************************************************" && \
	$(MAKE)); \
	if [ $$? -ne 0 ]; then \
		echo "!!!!!! ERROR !!!!!!" ; exit 1; \
	fi \
	done

clean_pzcl_kernel: $(PZCL_KERNEL_DIRS)
	@for d in $(PZCL_KERNEL_DIRS); \
	do \
	(cd $$d && \
	echo "" && \
	echo "*******************************************************" && \
	echo "* clean in ./$$d" && \
	echo "*******************************************************" && \
	$(MAKE) clean); \
	if [ $$? -ne 0 ]; then \
		echo "!!!!!! ERROR !!!!!!" ; exit 1; \
	fi \
	done


-include $(DEPS)


