# (i) Choose code
#para_type = serial
#para_type = omp
para_type = flat_mpi
#para_type = hybrid

# (ii) Compiler, compile options, etc.
FFTW_LOC  = /opt/klocal/fftw-3.3.6pl2-fujitsu
FFTW_INC  = -I$(FFTW_LOC)/include
ifeq ($(para_type),serial)
FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3 
CXX=time -p FCCpx
MACRO = -DPARALLELIZATION_TYPE=0
endif
ifeq ($(para_type),omp)
FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3_omp -lfftw3
CXX=time -p FCCpx -Kopenmp
MACRO = -DPARALLELIZATION_TYPE=1
endif
ifeq ($(para_type),flat_mpi)
FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3_mpi -lfftw3
CXX=time -p mpiFCCpx 
MACRO = -DPARALLELIZATION_TYPE=2
endif
ifeq ($(para_type),hybrid)
FFTW_LIB  = -L$(FFTW_LOC)/lib64 -lfftw3_mpi -lfftw3_omp -lfftw3
CXX=time -p mpiFCCpx -Kopenmp
MACRO = -DPARALLELIZATION_TYPE=3
endif
LDFLAGS   = $(FFTW_LIB)

#CXXFLAGS = -std=c++11 -O0 -g -NRtrap -Nquickdbg=subchk -Nquickdbg=heapchk -Nquickdbg=inf_detail
CXXFLAGS = -std=c++11 -Kfast -Ksimd=2 -Krestp=all
CXXFLAGS += $(MACRO) $(FFTW_INC)

# (iii) Job submission variables
QSUB = pjsub
ifeq ($(para_type),serial)
JOB_FILE_NAME = job.K.serial.sh
WORKDIR = /data/ra000008/namekata/PMMM/benchmark/fftw/serial
endif
ifeq ($(para_type),omp)
JOB_FILE_NAME = job.K.omp.sh
WORKDIR = /data/ra000008/namekata/PMMM/benchmark/fftw/omp
endif
ifeq ($(para_type),flat_mpi)
JOB_FILE_NAME = job.K.flat_mpi.sh
WORKDIR = /data/ra000008/namekata/PMMM/benchmark/fftw/flat_mpi
endif
ifeq ($(para_type),hybrid)
JOB_FILE_NAME = job.K.hybrid.sh
WORKDIR = /data/ra000008/namekata/PMMM/benchmark/fftw/hybrid
endif

# (iv) Executable
TARGET	= fft.out


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
