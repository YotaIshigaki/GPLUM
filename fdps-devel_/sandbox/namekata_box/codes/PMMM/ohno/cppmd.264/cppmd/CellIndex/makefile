
VPATH += ../Common 
INCDIRS += -I../Common

CXXFLAGS += -O3
#CXXFLAGS += --param max-inline-insns-single=4500
CXXFLAGS += -ffast-math -DUSE_TABLE_FOR_EWALD_REAL -g $(INCDIRS)
#CXXFLAGS += -Winline
CXXDEBUG = -g $(INCDIRS)
INCS += HalfShell.h CellIndex.h Common.h SpaceVector.h

testHalfShell : testHalfShell.o HalfShell.o CellIndex.o
	$(CXX) $(CXXFLAGS) -o $@ testHalfShell.o  HalfShell.o CellIndex.o

testHalfShell_debug : testHalfShell_debug.o HalfShell_debug.o CellIndex_debug.o
	$(CXX) $(CXXDEBUG) -o $@ testHalfShell_debug.o HalfShell_debug.o CellIndex_debug.o

testHalfShell_debug.o : testHalfShell.cc
	$(CXX) $(CXXDEBUG) -c testHalfShell.cc -o $@
HalfShell_debug.o : HalfShell.cc
	$(CXX) $(CXXDEBUG) -c $< -o $@
CellIndex_debug.o : CellIndex.cc
	$(CXX) $(CXXDEBUG) -c $< -o $@

testHalfShell.o : testHalfShell.cc $(INCS)

clean :
	rm -f *.o testHalfShell 
