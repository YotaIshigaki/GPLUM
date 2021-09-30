# (i) PATH setting
FFTW_LOC  = /opt/klocal/fftw-3.3.6pl2-fujitsu
FFTW_INC  = -I$(FFTW_LOC)/include
FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3
#FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3_omp -lfftw3
#FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3_mpi -lfftw3_omp -lfftw3

# (i) compiler, compile options
CXX = time -p mpiFCCpx
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all -Kopenmp $(FFTW_INC)
LDFLAGS = $(FFTW_LIB)

# (iv) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/PMMM/tests/fftw/

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
