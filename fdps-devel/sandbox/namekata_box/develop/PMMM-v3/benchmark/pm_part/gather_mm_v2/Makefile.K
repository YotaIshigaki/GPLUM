# Compiler, compile options, etc.
FFTW_LOC  = /opt/klocal/fftw-3.3.6pl2-fujitsu
FFTW_INC  = -I$(FFTW_LOC)/include
FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3_mpi -lfftw3 -lm

CXX=time -p mpiFCCpx 
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all $(FFTW_INC)
LDFLAGS = $(FFTW_LIB)

# (iii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
#WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pm_part/gather_mm_v2
#WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pm_part/gather_mm_v2/square_2d_no_offset
WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pm_part/gather_mm_v2/square_2d_offset_xy
#WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pm_part/gather_mm_v2/manhattan_2d_no_offset
#WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pm_part/gather_mm_v2/manhattan_2d_offset_xy
#WORKDIR = /data/ra000008/namekata/PMMM/benchmark/pm_part/gather_mm_v2/manhattan_3d

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
