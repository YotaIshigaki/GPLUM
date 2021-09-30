# (i) ffte 
FFTE_LOC = ./ffte-6.0
FFTE_INC = -I$(FFTE_LOC) -I$(FFTE_LOC)/mpi

# (i) Compiler, compile options, etc.
F77	= time -p mpifrtpx
FC	= time -p mpifrtpx
F77FLAGS = -Kfast -Ksimd=2 -Kopenmp $(FFTE_INC)
FCFLAGS = -Kfast -Ksimd=2 -Kopenmp $(FFTE_INC)
LDFLAGS = 

# (ii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/PMMM/benchmark/ffte

# (iii) Source files and executable
SRC_FFTE_SERIAL = $(FFTE_LOC)/dzfft2d.f \
		  $(FFTE_LOC)/dzfft3d.f \
		  $(FFTE_LOC)/factor.f  \
		  $(FFTE_LOC)/fft235.f  \
		  $(FFTE_LOC)/kernel.f  \
		  $(FFTE_LOC)/mfft235.f \
		  $(FFTE_LOC)/zdfft2d.f \
		  $(FFTE_LOC)/zdfft3d.f \
		  $(FFTE_LOC)/zfft1d.f  \
		  $(FFTE_LOC)/zfft2d.f  \
		  $(FFTE_LOC)/zfft3d.f
SRC_FFTE_MPI	= $(FFTE_LOC)/mpi/pdzfft2d.f \
		  $(FFTE_LOC)/mpi/pdzfft3d.f \
		  $(FFTE_LOC)/mpi/pfactor.f \
		  $(FFTE_LOC)/mpi/pzdfft2d.f \
		  $(FFTE_LOC)/mpi/pzdfft3d.f \
		  $(FFTE_LOC)/mpi/pzfft1d.f \
		  $(FFTE_LOC)/mpi/pzfft2d.f \
		  $(FFTE_LOC)/mpi/pzfft3d.f \
		  $(FFTE_LOC)/mpi/pzfft3dv.f
SRC_F77 	= $(SRC_FFTE_SERIAL) $(SRC_FFTE_MPI)
SRC_FC		= main.F90
OBJ_F77		= $(notdir $(SRC_F77:f=o))
OBJ_FC		= $(SRC_FC:F90=o) 
OBJ		= $(OBJ_F77) $(OBJ_FC)
VPATH	= $(FFTE_LOC) $(FFTE_LOC)/mpi
TARGET	= ffte.out

# Default rule
%.o : %.f
	$(F77) $(F77FLAGS) -c $<
%.o : %.F90
	$(FC) $(FCFLAGS) -c $<


all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
	mkdir -p $(WORKDIR)
	mv $@ $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

clean:
	rm -f *.o

distclean: clean
	rm -f $(TARGET)
	rm -f *.txt
