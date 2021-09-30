
VPATH += ../Common ../LJParameter ../SetPairList
INCDIRS += -I../Common -I../LJParameter -I../SetPairList

CXXFLAGS += -O2
CXXFLAGS += -Winline
#CXXFLAGS += --param max-inline-insns-single=4500
CXXFLAGS += -ffast-math -DUSE_TABLE_FOR_EWALD_REAL -g $(INCDIRS)
#CXXFLAGS += -Winline
INCS += ShortRangeInteraction.h ShortRangeInteractionSet.h ParticleInfo.h LJAmber94.h LJAmber.h LJParameter.h  Common.h SpaceVector.h

testShortRangeInteraction : testShortRangeInteraction.o EwaldRealInterpolation.o 
	$(CXX) $(CXXFLAGS) -o $@ testShortRangeInteraction.o EwaldRealInterpolation.o



testShortRangeInteraction.o : testShortRangeInteraction.cc $(INCS)

testPair: testPair.cc
	$(CXX) $(CXXFLAGS) -o testPair testPair.cc

clean :
	rm -f *.o testShortRangeInteraction testPair

callgrind$(TESTTYPE) : 
	valgrind --tool=callgrind --dump-instr=yes --trace-jump=yes ./testShortRangeInteraction -l 100 -t $(TESTTYPE)
