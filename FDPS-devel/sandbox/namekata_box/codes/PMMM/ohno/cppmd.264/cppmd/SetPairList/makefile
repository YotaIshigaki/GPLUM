
VPATH += ../Common ../LJParameter ../SetPairList
INCDIRS += -I../Common -I../LJParameter -I../SetPairList

CXXFLAGS += -O3 -ffast-math -DUSE_TABLE_FOR_EWALD_REAL -g -Winline $(INCDIRS)

INCS += ShortRangeInteraction.h ParticleInfo.h LJAmber94.h LJAmber.h LJParameter.h  Common.h SpaceVector.h

testPair: testPair.cc
	$(CXX) $(CXXFLAGS) -o testPair testPair.cc

