#CXX=mpiicpc
CXX=mpic++
VPATH += ../Common ../LJParameter ../SetPairList ../CellIndex
INCDIRS += -I../Common -I../LJParameter -I../SetPairList -I../CellIndex -I../Geometry -I../Timer -I../LongRangeInteraction -I../CalculationUnit -I../MPIParallel -I../ShortRangeInteraction -I../CovalentBond -I.. -I../convert -I../Dump -I../TestData -I/usr/local/boost/include/boost-1_41 -I../../../fftw2_gcc/include

#CXXFLAGS += -O3
CXXFLAGS += -O0
#CXXFLAGS += -ffast-math 
CXXFLAGS += -DUSE_TABLE_FOR_EWALD_REAL -g $(INCDIRS)
CXXFLAGS +=-DUSEMPI
CXXFLAGS +=-DPARTICLE_ISEND
#CXXFLAGS += -Winline
INCS += MPIParallel.h ParticleInfo.h CovalentBondInfo.h LJAmber94.h LJAmber.h LJParameter.h  SetPairList.h Common.h SpaceVector.h CubicCell.h

#FFTWLIB = -L/home/data11/ohno/md_code/fftw2_icc/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
FFTWLIB = -L../../../fftw2_gcc/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw

testMPIParallelLongRange : testMPIParallelLongRange.o MPIParallelLongRange.o  ../libcppmd.a
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm $(FFTWLIB) 

testMPIParallel : testMPIParallel.o ../libcppmd.a
	$(CXX) $(CXXFLAGS) -o $@ testMPIParallel.o ../libcppmd.a

testMove : testMove.o 
	$(CXX) $(CXXFLAGS) -o $@ testMove.o

testMPIParallel.o : testMPIParallel.cc $(INCS)
testMove.o : testMove.cc $(INCS)

clean:
	rm -f *.o testMPIParallel testMPIParallelLongRange

testsamplecell: testsamplecell.cc samplecell.h samplecell.o
	$(CXX) $(CXXFLAGS) -o testsamplecell testsamplecell.cc samplecell.o

../libcppmd.a :
	cd ../..; make



