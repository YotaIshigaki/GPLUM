CXX=g++
CXXFLAGS+=-O3 -g
#CXXFLAGS+=-DUSE_PAPI
#PAPI_LIB = /usr/local/lib/libpapi.a

testTimer : testTimer.cc Timer.h Timer.o
	$(CXX) $(CXXFLAGS) -o $@ $< Timer.o $(PAPI_LIB)

Timer.o : Timer.cc Timer.h
	$(CXX) $(CXXFLAGS) -c $<

clean :
	rm -f testTimer *.o
