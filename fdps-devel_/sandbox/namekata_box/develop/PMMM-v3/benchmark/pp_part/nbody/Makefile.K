# (i) PATH setting
PS_PATH = ../FDPS-5.0d/src/
INC = -I$(PS_PATH)

# (ii) Compiler and compile options
CXX = time -p mpiFCCpx
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# (iii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
#WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pp_part/1M
#WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pp_part/100k
#WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pp_part/10k
WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pp_part/1k

# (iv) Executable
TARGET	= nbody.out

all:$(TARGET)

$(TARGET): nbody.cpp $(OBJS)
	$(CXX) $(INC) $(CXXFLAGS) -o $@ $^ $(CLIBS)
	mkdir -p $(WORKDIR)
	mv $(TARGET) $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

clean:
	rm -f *.o *~

distclean: clean
	rm -f nbody.out
	rm -rf result
