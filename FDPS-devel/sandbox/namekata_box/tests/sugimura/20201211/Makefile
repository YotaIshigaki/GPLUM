PS_PATH = -I ../../../../../src/

CXX = time -p mpiicpc
CXXFLAGS = -std=c++11 -O3 -ipo -xCORE-AVX512 -no-prec-div
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -qopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

QSUB = pjsub
JOB_FILE_NAME = job.ito-subsys-a.sh
WORKDIR = ./work

OBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HDRS = $(wildcard *.h)
PROGRAM = test.out

.PHONY:	clean all

all:	$(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(PROGRAM) $(PS_PATH)
	mkdir -p $(WORKDIR)
	mv $(PROGRAM) $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

%.o:	%.cpp $(HDRS)
	$(CXX) -c $< $(CXXFLAGS) $(PS_PATH)

clean:
	rm -f $(PROGRAM) $(OBJS)
