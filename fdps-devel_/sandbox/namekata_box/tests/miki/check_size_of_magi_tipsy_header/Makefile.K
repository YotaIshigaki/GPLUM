# Compiler, compile options, etc.
FC=time -p frtpx
CXX=time -p FCCpx 
FCFLAGS = -X 03 -Kfast -Ksimd=2
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
LDFLAGS = 

# (iii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/tests/check_size_of_magi_tipsy_header


all: 
	$(CXX) $(CXXFLAGS) -o cpp.out main.cpp $(LDFLAGS)
	$(FC)  $(FCFLAGS)  -o ftn.out main.F90 $(LDFLAGS)
	mkdir -p $(WORKDIR)
	mv *.out $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

clean:
	rm -f *.o

distclean: clean
	rm -f $(TARGET)
	rm -f *.txt
