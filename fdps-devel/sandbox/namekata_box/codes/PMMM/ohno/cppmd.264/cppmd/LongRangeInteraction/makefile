
VPATH += ../Common
INCDIRS += -I../Common -I.. -I../..
INCDIRS += -I../convert
VPATH += ../Geometry
INCDIRS += -I../Geometry


CXXFLAGS += -g -O2 $(INCDIRS) -UNDEBUG
CXXFLAGS += -DHAVE_CONFIG_H
CXXFLAGS += $(FFTFLAGS)
CXXFLAGS += -I/home/data11/ohno/md_code/fftw2/include -L/home/data11/ohno/md_code/fftw2/lib
CXXFLAGS += -DLONGRANGE_DT
#CXXFLAGS += -DHAVE_CONFIG_H -DNOT_USE_MPI
#CXXFLAGS += -m32
LDFLAGS += -lrt

INCS += Common.h SpaceVector.h ParticleInfo.h LongRangeInteraction.h RealSpaceLongRangeInteraction.h
OBJS = testLongRangeInteraction.o RealSpaceLongRangeInteraction.o
PMEOBJS = PMELongRangeInteraction.o \
          PMEMethod.o FFT3D.o PMEArray.o ParallelPMEMethod.o EwaldBase.o \
          Geometry.o
EWALDOBJS = EwaldLongRangeInteraction.o \
          EwaldMethod.o ParallelEwald.o EwaldBase.o \
          Geometry.o
MGOBJS = MultigridLongRangeInteraction.o \
         LGM.o MGPoisson.o FFT3D.o Geometry.o

all : longrange testLongRangeInteraction \
      testPMELongRangeInteraction testEwaldLongRangeInteraction testPMEEwald \
      testMultigridLongRangeInteraction testAnalyticalGrid testMultigrid

longrange : longrange.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lfftw -lm $(LDFLAGS)

testLongRangeInteraction : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lfftw -lm $(LDFLAGS)

#testPMELongRangeInteraction : FFTFLAGS = -DHAVE_CONFIG_H -DNOT_USE_MPI
testPMELongRangeInteraction : FFTFLAGS = -DUSE_FFTW2
testPMELongRangeInteraction : testPMELongRangeInteraction.o $(PMEOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -lm $(LDFLAGS)
testghostPMELongRangeInteraction : FFTFLAGS = -DUSE_FFTW2
testghostPMELongRangeInteraction : testghostPMELongRangeInteraction.o $(PMEOBJS) $(EWALDOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -lm $(LDFLAGS)

#testPMEEwald : FFTFLAGS = -DHAVE_CONFIG_H -DNOT_USE_MPI
testPMEEwald : FFTFLAGS = -DUSE_FFTW2
testPMEEwald : testPMEEwald.o $(PMEOBJS) $(EWALDOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lrfftw -lfftw -lm $(LDFLAGS)

testEwaldLongRangeInteraction : testEwaldLongRangeInteraction.o $(EWALDOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm $(LDFLAGS)

testMultigridLongRangeInteraction : testMultigridLongRangeInteraction.o $(MGOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lrfftw -lfftw -lm $(LDFLAGS)

testMultigrid : testMultigrid.o $(MGOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lrfftw -lfftw -lm $(LDFLAGS)

check : all
	./testLongRangeInteraction && ./longrange \
        && ./testPMEEwald \
	&& ./testPMELongRangeInteraction && ./testEwaldLongRangeInteraction \
	&& ./testMultigridLongRangeInteraction

mpi: CXX=mpicxx -DUSE_MPI -DUSE_FFTW2
mpi: cleanobj testMultigridMPI testPMEEwaldMPI
testMultigridMPI : testMultigridMPI.o $(MGOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lrfftw -lfftw -lrfftw_mpi -lfftw_mpi -lm $(LDFLAGS)
testPMEEwaldMPI : testPMEEwaldMPI.o $(PMEOBJS) $(EWALDOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lrfftw -lfftw -lrfftw_mpi -lfftw_mpi -lm $(LDFLAGS)


cleanobj :
	-rm -f *.o

clean :
	-rm -f longrange.o longrange $(OBJS) testLongRangeInteraction
	-rm -f $(PMEOBJS) testPMELongRangeInteraction{.o,}
	-rm -f $(EWALDOBJS) testEwaldLongRangeInteraction{.o,}
	-rm -f testPMEEwald{.o,}
	-rm -f testAnalyticalGrid
	-rm -f $(MGOBJS) testMultigridLongRangeInteraction{.o,}
	-rm -f testMultigrid{.o,}
	-rm -f testPMEEwaldMPI{.o,}
	-rm -f testMultigridMPI{.o,}
