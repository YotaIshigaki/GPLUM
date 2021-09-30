# (i) PATH setting
PS_PATH = ./src/
INC = -I$(PS_PATH)

# (i) Compiler, compile options, etc.
#CXX=time -p FCCpx
CXX=time -p mpiFCCpx 
#CXXFLAGS = -std=c++11 -O0 -g -NRtrap -Nquickdbg=subchk -Nquickdbg=heapchk -Nquickdbg=inf_detail
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DREAD_DATA_WITH_BYTESWAP
LDFLAGS	= 

# (ii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
#WORKDIR = /data/ra000008/namekata/GB2020/cluster_identification/30kyr
#WORKDIR = /data/ra000008/namekata/GB2020/cluster_identification/20kyr
WORKDIR = /data/ra000008/namekata/GB2020/cluster_identification/10kyr

# (iii) Data preprocess commands
MKDATA = ./scripts/data_partition.py
ORG_DATA_DIR = /data/ra000008/namekata/GB2020/cluster_identification/input_data_org
#ORG_DATA_FILE_NAME = snapshot_0067
#ORG_DATA_FILE_NAME = snapshot_0155
ORG_DATA_FILE_NAME = snapshot_1001
ORG_DATA_FILE = $(ORG_DATA_DIR)/$(ORG_DATA_FILE_NAME)
DATADIR = $(WORKDIR)/input_data

# (iv) Executable
TARGET	= a.out

all: $(TARGET)

$(TARGET): main.o
	$(CXX) $(INC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	mkdir -p $(WORKDIR)
	mv $@ $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	$(MKDATA) $(JOB_FILE_NAME) $(ORG_DATA_FILE) $(DATADIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

main.o: main.cpp 
	$(CXX) $(INC) $(CXXFLAGS) -c $< 

clean:
	rm -f *.o

distclean: clean
	rm -f $(TARGET)
	rm -f *.txt
