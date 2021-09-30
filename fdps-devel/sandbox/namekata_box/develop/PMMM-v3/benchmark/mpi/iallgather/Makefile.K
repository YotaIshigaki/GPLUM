# Compiler, compile options, etc.
CXX=time -p mpiFCCpx 
#CXXFLAGS = -std=c++11 -O0 -g -NRtrap -Nquickdbg=subchk -Nquickdbg=heapchk -Nquickdbg=inf_detail
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
LDFLAGS = 

# (iii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/PMMM/benchmark/mpi/iallgather

# (iv) Executable
TARGET	= iallgather.out


all: $(TARGET)

$(TARGET): main.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	mkdir -p $(WORKDIR)
	mv $@ $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

main.o: main.cpp 
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -f *.o

distclean: clean
	rm -f $(TARGET)
	rm -f *.txt
