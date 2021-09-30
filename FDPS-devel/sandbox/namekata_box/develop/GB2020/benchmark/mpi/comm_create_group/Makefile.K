# (i) Compiler, compile options, etc.
CXX=time -p mpiFCCpx 
#CXXFLAGS = -std=c++11 -O0 -g -NRtrap -Nquickdbg=subchk -Nquickdbg=heapchk -Nquickdbg=inf_detail
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
LDFLAGS	= 

# (ii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/GB2020/benchmark/mpi/comm_create_group

# (iv) Executable
TARGET	= a.out

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
