LDFLAGS   = 
CXX = time -p FCCpx
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = a.out

.PHONY:	clean all

all: $(CPPOBJS) $(CPPHDRS)
	$(CXX) $(CXXFLAGS) $(CPPOBJS) -o $(PROGRAM) $(LDFLAGS)

%.o: %.cpp $(CPPHDRS)
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -f $(CPPOBJS)

distclean: clean
	rm -f $(PROGRAM)
