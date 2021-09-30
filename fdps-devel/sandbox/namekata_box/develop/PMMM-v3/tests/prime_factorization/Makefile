# (i) compiler, compile options
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O3
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
