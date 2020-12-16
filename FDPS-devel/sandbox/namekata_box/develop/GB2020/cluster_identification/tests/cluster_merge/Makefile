# (i) FDPS
FDPS_LOC = ../../
FDPS_INC = -I$(FDPS_LOC)/src

# (i) compiler, compile options
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O0 $(FDPS_INC)
LDFLAGS = 

TARGET = a.out

all: main.o
	$(CXX) $(CXXFLAGS) $< -o $(TARGET) $(LDFLAGS)

main.o: main.cpp
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -f *.o

distclean: clean
	rm -f $(TARGET)
