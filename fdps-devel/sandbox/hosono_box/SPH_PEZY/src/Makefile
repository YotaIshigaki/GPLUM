PZCL_KERNEL_DIRS=kernel.sc32

INC = -I ../FDPS/src
#CC = clang++
CC = mpic++
CFLAGS = -O3 #-ffast-math -funroll-loops 
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith
CFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
CFLAGS += -DFAST_COMM_FOR_2ND_EXCHANGE
#CFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

CFLAGS += -DENABLE_PEZY

PZSDK_PATH = /opt/pzsdk.ver3.0
INC += -I$(PZSDK_PATH)/inc
PZSDK_LIB_PATH = $(PZSDK_PATH)/lib

PZCL_STATIC_LIB =1
ifeq ($(PZCL_STATIC_LIB),1)
PZCL_LIB = $(PZSDK_LIB_PATH)/libpzcl.a
else
PZCL_LIB = -lpzcl
endif

LIBS += $(PZCL_LIB) -L$(PZSDK_LIB_PATH) -lpthread -ldl -lstdc++ -lm

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = SPzyH.out

.PHONY:	clean all

all:	make_pzcl_kernel $(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(INC)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Building $< ..."
	@$(CC) -c $< $(CFLAGS) $(WARNINGS) $(INC)
	@echo "[$< OK]"

clean:
	-rm *.out *.o *.log

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

