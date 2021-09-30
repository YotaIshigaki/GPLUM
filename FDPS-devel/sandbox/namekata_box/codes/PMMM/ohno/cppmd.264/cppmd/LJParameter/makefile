
VPATH += ../Common ../LJParameter 
INCDIRS += -I../Common -I../LJParameter 

CXXFLAGS += -O3 -ffast-math -g $(INCDIRS)
#CXXFLAGS += -Winline
TARGET=LJAmber94
INCS += $(TARGET).h LJAmber.h LJParameter.h  Common.h SpaceVector.h

test$(TARGET) : test$(TARGET).o LJParameter.o
	$(CXX) $(CXXFLAGS) -o $@ test$(TARGET).o  LJParameter.o



test$(TARGET).o : test$(TARGET).cc $(INCS)


