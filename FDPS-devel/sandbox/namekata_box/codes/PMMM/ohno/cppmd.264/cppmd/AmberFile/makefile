
VPATH += ../Common 
INCDIRS += -I../Common

CXXFLAGS += -g $(INCDIRS) -UNDEBUG
#CXXFLAGS += -m32
INCS += Common.h SpaceVector.h CovalentBondInfo.h AmberFile.h AmberPrmtop.h AmberRestrt.h ReadAmber.h ReadAmber-inl.h debug.hh

testAmberFile : testAmberFile.o AmberFile.o AmberPrmtop.o AmberRestrt.o ReadAmber.o debug.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean :
	rm -f *.o testAmberFile
