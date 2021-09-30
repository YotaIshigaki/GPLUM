# Compiler, compile options, etc.
CXX=time -p FCCpx 
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
LDFLAGS = 

# (iii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/tests/io_binary_file


all: 
	# Make data.dat
	g++ -o write.out write.cpp 
	./write.out
	# Perform read.out
	$(CXX) $(CXXFLAGS) -o read.out read.cpp $(LDFLAGS)
	mkdir -p $(WORKDIR)
	mv read.out $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cp data.dat $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

clean:
	rm -f *.o

distclean: clean
	rm -f $(TARGET)
	rm -f *.txt
