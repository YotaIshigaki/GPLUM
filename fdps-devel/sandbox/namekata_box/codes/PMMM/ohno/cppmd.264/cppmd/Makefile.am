## Copyright (c) 2008, 2009, 2010, 2011, 2012 RIKEN. All Rights Reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##    1. Redistributions of source code must retain the above copyright
##    notice, this list of conditions and the following disclaimer.
##    2. Redistributions in binary form must reproduce the above
##    copyright notice, this list of conditions and the following
##    disclaimer in the documentation and/or other materials provided
##    with the distribution.
##    3. The name of the author may not be used to endorse or promote
##    products derived from this software without specific prior written
##    permission.
##
##
## THIS SOFTWARE IS PROVIDED BY RIKEN ``AS IS'' AND ANY EXPRESS OR
## IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
## HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
## STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
## IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.

DISTCLEANFILES = cppmd-int.h
CLEANFILES = ewald_table.mod shift_function.mod

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/cppmd/AmberFile -I$(top_srcdir)/cppmd/CellIndex -I$(top_srcdir)/cppmd/Common -I$(top_srcdir)/cppmd/CovalentBond -I$(top_srcdir)/cppmd/LJParameter -I$(top_srcdir)/cppmd/CalculationUnit -I$(top_srcdir)/cppmd/LongRangeInteraction -I$(top_srcdir)/cppmd/Geometry -I$(top_srcdir)/cppmd/MPIParallel -I$(top_srcdir)/cppmd/SetPairList -I$(top_srcdir)/cppmd/ShortRangeInteraction -I$(top_srcdir)/cppmd/PairList -I$(top_srcdir)/cppmd/TestData -I$(top_srcdir)/cppmd/Timer -I$(top_srcdir)/cppmd/Dump -I$(top_srcdir)/cppmd/HDF -I$(top_srcdir)/cppmd/PMMM $(BOOST_CPPFLAGS) $(HDFINC)

LONGRANGESOURCES = \
  LongRangeInteraction/LongRangeParameter.h \
  LongRangeInteraction/EwaldBase.cc \
  LongRangeInteraction/EwaldBase.h \
  LongRangeInteraction/NullStream.h \
  LongRangeInteraction/UseConfig.h \
  LongRangeInteraction/UseMPI.h

MR3EXAFMMSOURCES = \
  LongRangeInteraction/MR3_FMM.cc \
  LongRangeInteraction/MR3_FMM.h
EXAFMMSOURCES = \
  EXAFMM/cpuCartesianLaplace.h \
  EXAFMM/cpuEvaluator.h \
  EXAFMM/cpuP2P.h \
  EXAFMM/cpuSphericalLaplace.h \
  EXAFMM/bottomup.h \
  EXAFMM/evaluator.h \
  EXAFMM/kernel.h \
  EXAFMM/logger.h \
  EXAFMM/mympi.h \
  EXAFMM/parallelfmm.h \
  EXAFMM/partition.h \
  EXAFMM/serialfmm.h \
  EXAFMM/sort.h \
  EXAFMM/tree.h \
  EXAFMM/types.h \
  EXAFMM/vec.h
EXAFMM2SOURCES = \
  EXAFMM2/evaluator.h \
  EXAFMM2/kernel.h \
  EXAFMM2/kernelcore.h \
  EXAFMM2/parallelfmm.h \
  EXAFMM2/serialfmm.h \
  EXAFMM2/types.h

PMMMSOURCES = \
  PMMM/PMMMInteraction.h \
  PMMM/pmmm_fmm.h \
  PMMM/M2L_convolution_PBC.h \
  PMMM/GreenFunction_PBC.h \
  PMMM/cutoff.h \
  MPIParallel/MPIParallelLongRangePMMM.h \
  CalculationUnit/PMMMPreparator.h


SIMPLEFFTSOURCES = \
  LongRangeInteraction/FFT.cc \
  LongRangeInteraction/FFT.h 

OLDFFTSOURCES = \
  LongRangeInteraction/FFT3D.cc \
  LongRangeInteraction/FFT3D.h \
  LongRangeInteraction/FFT3Dffte.h \
  LongRangeInteraction/FFT3Dfftw2.h \
  LongRangeInteraction/FFT3Dfftw3.h \
  LongRangeInteraction/FFT3Dmkl.h \
  LongRangeInteraction/MPIFFT3D.h \
  LongRangeInteraction/MPIFFT3Dffte.h \
  LongRangeInteraction/MPIFFT3Dfftw2.h \
  LongRangeInteraction/MPIFFT3Dmkl.h 

PMESOURCES = \
  LongRangeInteraction/Array3D.h \
  LongRangeInteraction/PMEArray.cc \
  LongRangeInteraction/PMEArray.h \
  LongRangeInteraction/PME.cc \
  LongRangeInteraction/PME.h

OLDPMESOURCES = \
  LongRangeInteraction/Array3D.h \
  LongRangeInteraction/ParallelPMEMethod.cc \
  LongRangeInteraction/PMEArray.cc \
  LongRangeInteraction/PMEArray.h \
  LongRangeInteraction/PMEInterface.h \
  LongRangeInteraction/PMEInterfaceFwd.h \
  LongRangeInteraction/PMEInterfaceImpl.h \
  LongRangeInteraction/PMELongRangeInteraction.cc \
  LongRangeInteraction/PMELongRangeInteraction.h \
  LongRangeInteraction/PMEMethod.cc \
  LongRangeInteraction/PMEMethod.h

EWALDSOURCES = \
  LongRangeInteraction/EwaldLongRangeInteraction.cc \
  LongRangeInteraction/EwaldMethod.cc \
  LongRangeInteraction/RealSpaceLongRangeInteraction.cc \
  LongRangeInteraction/ParallelEwald.cc
EWALDHEADERS = \
  LongRangeInteraction/BSpline.h \
  LongRangeInteraction/Lagrange.h \
  LongRangeInteraction/EwaldInterface.h \
  LongRangeInteraction/EwaldInterfaceBase.h \
  LongRangeInteraction/EwaldInterfaceFwd.h \
  LongRangeInteraction/EwaldInterfaceImpl.h \
  LongRangeInteraction/EwaldLongRangeInteraction.h \
  LongRangeInteraction/EwaldMethod.h \
  LongRangeInteraction/RealSpaceLongRangeInteraction.h \
  LongRangeInteraction/SurfaceDipole.h

STMEWALDSOURCES = \
  LongRangeInteraction/StMELongRangeInteraction.cc \
  LongRangeInteraction/StMELongRangeInteraction.h

MULTIGRIDSOURCES = \
  LongRangeInteraction/MGPoisson.cc \
  LongRangeInteraction/MGPoisson.h \
  LongRangeInteraction/MultigridFwd.h \
  LongRangeInteraction/MultigridLongRangeInteraction.cc \
  LongRangeInteraction/MultigridLongRangeInteraction.h \
  LongRangeInteraction/mgp_layout_struct.h \
  LongRangeInteraction/ChargeAssignmentMethod.h \
  LongRangeInteraction/LGM.cc \
  LongRangeInteraction/LGM.h \
  LongRangeInteraction/MGFFT3D.h \
  LongRangeInteraction/MGFFT3DImpl.h

GFMMSOURCES = \
  FMM/CommMEC.cpp \
  FMM/FMMSolver.cpp \
  FMM/FMMSolver.h \
  FMM/cell_routines.cpp \
  FMM/cell.cpp \
  FMM/cell.h \
  FMM/coefficients.cpp \
  FMM/coefficients.h \
  FMM/debug_routines.cpp \
  FMM/defs.h \
  FMM/downward_pass.cpp \
  FMM/expansion_routines.cpp \
  FMM/io_routines.cpp \
  FMM/math_routines.cpp \
  FMM/upward_pass.cpp \
  LongRangeInteraction/FMMLongRangeInteraction.h

HDFSOURCES = \
  HDF/HDFSnapshotFile.cc \
  HDF/HDFSnapshotFile.h \
  HDF/HDFSnapshotData.cc \
  HDF/HDFSnapshotData.h \
  Dump/Bin2HDF.cc \
  Dump/BIn2HDF.h \
  Dump/HDF2RST.cc \
  Dump/HDF2RST.h \
  Dump/HDF2Bin.cc \
  Dump/HDF2Bin.h

bin_PROGRAMS = cppmd Bin2RST
cppmd_SOURCES = main.cc Cppmd.cc Cppmd.h Config.cc Config.h
cppmd_LDFLAGS = $(BOOST_LDFLAGS) $(HDF_LDFLAGS)
cppmd_LDADD = libcppmd.a $(HDF_LDADD)
Bin2RST_SOURCES = Dump/mainBin2RST.cc
Bin2RST_LDFLAGS = $(BOOST_LDFLAGS) 
Bin2RST_LDADD = libcppmd.a

Bin2HDF_SOURCES = Dump/mainBin2HDF.cc
Bin2HDF_LDFLAGS = $(BOOST_LDFLAGS) $(HDF_LDFLAGS)
Bin2HDF_LDADD = libcppmd.a $(HDF_LDADD)

HDF2RST_SOURCES = Dump/mainHDF2RST.cc
HDF2RST_LDFLAGS = $(BOOST_LDFLAGS) $(HDF_LDFLAGS)
HDF2RST_LDADD = libcppmd.a $(HDF_LDADD)

HDF2Bin_SOURCES = Dump/mainHDF2Bin.cc
HDF2Bin_LDFLAGS = $(BOOST_LDFLAGS) $(HDF_LDFLAGS)
HDF2Bin_LDADD = libcppmd.a $(HDF_LDADD)


lib_LIBRARIES = libcppmd.a

libcppmd_a_CXXFLAGS = -DPARTICLE_ISEND
libcppmd_a_SOURCES = \
  simulator.cc simulator.h \
  getline.cc getline.h \
  cppmd-int.h \
  \
  AmberFile/AmberFile.cc AmberFile/AmberFile.h \
  AmberFile/AmberPrmtop.cc AmberFile/AmberPrmtop.h \
  AmberFile/AmberRestrt.cc AmberFile/AmberRestrt.h \
  AmberFile/ReadAmber.cc AmberFile/ReadAmber.h \
  AmberFile/ReadAmber-inl.h \
  \
  CellIndex/CellIndex.cc CellIndex/CellIndex.h \
  CellIndex/CubicCell.cc CellIndex/CubicCell.h \
  CellIndex/HalfShell.cc CellIndex/HalfShell.h \
  CellIndex/CellMethod.h \
  CellIndex/CellMethodCore.h \
  CellIndex/CellMethodFwd.h \
  CellIndex/CellMethodImpl.h \
  CellIndex/HalfShellCore.h \
  CellIndex/HalfShellCoreImpl.h \
  CellIndex/SmallBallCore.h \
  CellIndex/SmallBallCoreImpl.h \
  \
  Common/Common.h \
  Common/CovalentBondInfo.h \
  Common/ErrorPos.h \
  Common/OperationSelector.h \
  Common/ParallelCoordinator.h \
  Common/ParticleInfo.h \
  Common/SpaceVector.h \
  Common/UnitParameter.h \
  Common/Log.h \
  \
  CovalentBond/CovalentBondInfo.cc \
  CovalentBond/CovalentBond.h \
  CovalentBond/CBInterface.h \
  CovalentBond/CBInterfaceFwd.h \
  CovalentBond/CBInterfaceImpl.h \
  CovalentBond/CBUnordered.h \
  CovalentBond/CBObjects.h \
  CovalentBond/CBObjects.cc \
  \
  LJParameter/LJAmber.h \
  LJParameter/LJAmber94.h \
  LJParameter/LJAr.h \
  LJParameter/LJParameter.cc LJParameter/LJParameter.h \
  \
  LongRangeInteraction/DummyLongRangeInteraction.cc \
  LongRangeInteraction/DummyLongRangeInteraction.h \
  LongRangeInteraction/LongRangeInteraction.h \
  \
  Geometry/Geometry.cc Geometry/Geometry.h \
  \
  MPIParallel/MPIParallel.cc MPIParallel/MPIParallel.h \
  MPIParallel/MPIParallelLongRange.cc MPIParallel/MPIParallelLongRange.h \
  \
  SetPairList/SetPairList.cc SetPairList/SetPairList.h \
  \
  ShortRangeInteraction/ShortRangeInteraction.cc \
  ShortRangeInteraction/ShortRangeInteraction.h \
  ShortRangeInteraction/ShortRangeInteractionSet.cc \
  ShortRangeInteraction/ShortRangeInteractionSet.h \
  ShortRangeInteraction/ShortRangeInteraction_noinline.h \
  \
  PairList/PairList.cc \
  PairList/PairList.h \
  PairList/PairListInteraction.cc \
  PairList/PairListInteraction.h \
  PairList/PairListZeroDipole.cc \
  PairList/PairListZeroDipole.h \
  \
  TestData/NoInputSystem.h \
  TestData/ArFCC.cc TestData/ArFCC.h \
  TestData/NaClFCC.cc TestData/NaClFCC.h \
  TestData/WaterLattice.cc TestData/WaterLattice.h \
  TestData/TestParticle.cc TestData/TestParticle.h \
  \
  Timer/Timer.cc Timer/Timer.h\
  \
  Dump/Dump.cc Dump/Dump.h\
  Dump/Converter.cc Dump/Converter.h\
  Dump/Bin2RST.cc Dump/Bin2RST.h\
  \
  CalculationUnit/CalculationUnit.cc CalculationUnit/CalculationUnit.h \
  CalculationUnit/CalcPreparator.cc CalculationUnit/CalcPreparator.h \
  CalcForce.cc CalcForce.h \
  Integrator.cc Integrator.h \
  Settle.cc Settle.h \
  Log.cc

if ENABLE_MR3EXAFMM
libcppmd_a_SOURCES += $(LONGRANGESOURCES) $(EWALDHEADERS)
libcppmd_a_SOURCES += $(MR3EXAFMMSOURCES) $(EXAFMMSOURCES) $(EXAFMM2SOURCES)
endif
if ENABLE_PMMM
libcppmd_a_SOURCES += $(LONGRANGESOURCES) $(EWALDHEADERS)
libcppmd_a_SOURCES += $(PMMMSOURCES)
endif
if ENABLE_PME
libcppmd_a_SOURCES += $(LONGRANGESOURCES) $(EWALDHEADERS)
libcppmd_a_SOURCES += $(PMESOURCES)
if ENABLE_SIMPLE_FFT
libcppmd_a_SOURCES += $(SIMPLEFFTSOURCES)
else
libcppmd_a_SOURCES += $(OLDFFTSOURCES)
endif
endif
if ENABLE_OLDPME
libcppmd_a_SOURCES += $(LONGRANGESOURCES) $(EWALDHEADERS)
libcppmd_a_SOURCES += $(OLDPMESOURCES)
libcppmd_a_SOURCES += $(OLDFFTSOURCES)
endif
if ENABLE_EWALD
libcppmd_a_SOURCES += $(LONGRANGESOURCES)
libcppmd_a_SOURCES += $(EWALDSOURCES) $(EWALDHEADERS)
endif
#if ENABLE_STMEWALD
#libcppmd_a_SOURCES += $(LONGRANGESOURCES) $(EWALDHEADERS)
#libcppmd_a_SOURCES += $(STMEWALDSOURCES)
##libcppmd_a_SOURCES += $(MULTIGRIDSOURCES)
#endif
#if ENABLE_GFMM
#libcppmd_a_SOURCES += $(LONGRANGESOURCES) $(EWALDHEADERS)
#libcppmd_a_SOURCES += $(GFMMSOURCES)
#endif
if USE_HDF
bin_PROGRAMS += Bin2HDF HDF2RST HDF2Bin
libcppmd_a_SOURCES += $(HDFSOURCES)
endif


if ENABLE_TABLE_FOR_EWALD_REAL
libcppmd_a_SOURCES += ShortRangeInteraction/EwaldRealInterpolation.cc \
		      ShortRangeInteraction/EwaldRealInterpolation.h
endif

if ENABLE_FORTRAN_KERNEL
libcppmd_a_SOURCES += PairList/PairListInteraction_f.f90
endif

EXTRA_DIST = \
  AmberFile/debug.cc AmberFile/debug.h \
  \
  CellIndex/makefile CellIndex/testHalfShell.cc \
  \
  Common/testRange.cc \
  \
  LJParameter/makefile LJParameter/testLJAmber94.cc \
  LJParameter/testLJAr.cc \
  \
  MPIParallel/makefile MPIParallel/old_testMPIParallel.cc \
  MPIParallel/samplecell.cc MPIParallel/samplecell.h \
  MPIParallel/testMPIParallel.cc MPIParallel/testMove.cc \
  MPIParallel/testsamplecell.cc \
  \
  SetPairList/concept_SetPairList.txt SetPairList/makefile \
  SetPairList/memo.txt SetPairList/testPair.cc \
  \
  ShortRangeInteraction/EwaldRealInterpolation.cpp.200r.mach \
  ShortRangeInteraction/makefile \
  ShortRangeInteraction/testShortRangeInteraction.cc \
  \
  TestData/ar20.pdb TestData/makefile TestData/testArFCC.cc \
  \
  Doxyfile memo.txt todo \
  \
  EXAFMM2/Makefile \
  EXAFMM2/main.cxx \
  \
  $(STMEWALDSOURCES) \
  $(MULTIGRIDSOURCES) \
  $(GFMMSOURCES)

## vim:se ft=automake:
