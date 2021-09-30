
VPATH += ../Common 
INCDIRS += -I../Common

CXXFLAGS += -O3
#CXXFLAGS += --param max-inline-insns-single=4500
CXXFLAGS += -ffast-math -DUSE_TABLE_FOR_EWALD_REAL -g $(INCDIRS)
#CXXFLAGS += -Winline
CXXDEBUG = -g $(INCDIRS)
INCS += NoInputSystem.h ArFCC.h NaClFCC.h Common.h SpaceVector.h

all : testArFCC testNaClFCC

testArFCC : testArFCC.o ArFCC.o 
	$(CXX) $(CXXFLAGS) -o $@ testArFCC.o  ArFCC.o 

testNaClFCC : testNaClFCC.o NaClFCC.o
	$(CXX) $(CXXFLAGS) -o $@ testNaClFCC.o NaClFCC.o

clean :
	rm -f *.o testArFCC testNaClFCC
