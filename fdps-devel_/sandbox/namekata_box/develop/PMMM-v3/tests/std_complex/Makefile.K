# (i) compiler, compile options
CXX = FCCpx
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
LDFLAGS = 

# (iv) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/PMMM/tests/std_complex/

TARGET = a.out

all: main.o
	$(CXX) $(CXXFLAGS) $< -o $(TARGET) $(LDFLAGS)
	mkdir -p $(WORKDIR)
	mv $(TARGET) $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

main.o: main.cpp 
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -f *.o

distclean: clean
	rm -f $(TARGET)
