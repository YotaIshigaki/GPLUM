FDPS_PATH?=$(HOME)/FDPS
INC+= -I$(FDPS_PATH)/src

CXX=g++-9
CXXFLAGS= -O2 $(INC)
#CXXFLAGS+= -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CXXFLAGS+= -DPARTICLE_SIMULATOR_MPI_PARALLEL

SRC=main.cpp
HED=

all:	argon.out
argon.out:	$(SRC) Makefile
	$(CXX) $(CXXFLAGS) $(SRC) -o argon.out

test:	argon.out
	./argon.out 10

clean:
	rm argon.out
