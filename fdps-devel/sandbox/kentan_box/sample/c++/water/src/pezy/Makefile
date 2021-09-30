#======================================================================
#   Numerical Libraries and Compilers
#======================================================================

FDPS_LOC = $(HOME)/fdps_svn
FDPS_INC = -I$(FDPS_LOC)/src -I../../inc

CXXFLAGS_COMMON = -std=c++11 -O3 -ffast-math -funroll-loops
#CXXFLAGS_COMMON = -std=c++11 -xHost -ip

CXXFLAGS_COMMON += $(FDPS_INC)
CXXFLAGS_COMMON += -DENABLE_PEZY

PZCL_KERNEL_DIRS=pezy/kernel.sc32

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

# [1] Serial
CXX = g++
CXXFLAGS = $(CXXFLAGS_COMMON)
#CXXFLAGS+= -DPARTICLE_SIMULATOR_DEBUG_PRINT
# [2] OpenMP
#CXX = g++
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# [3] MPI
#CXX = mpic++
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -L$(MPI_PATH)/lib
# [4] MPI + OpenMP
CXX = mpic++
CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

CXXFLAGS += -DENABLE_REUSE
#----------------------------------------------------------------------
#   Source and Object files
#----------------------------------------------------------------------
SRC   = main.cpp
HED   = user_defined_class.h kernel.hpp

all:	pezy

pezy:	make_pzcl_kernel $(OBJ_DIR) Makefile $(SRC) $(HED)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SRC) $(LIBS) $(LIB_DIR) $(INC_DIR) -o pezy.out

clean:	clean_pzcl_kernel
	rm -f *.out *.o

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

