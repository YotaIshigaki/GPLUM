PS_PATH_FDPS = ../FDPS-5.0g/src/
PS_PATH_PIKG = ../PIKG/inc/
INC = -I$(PS_PATH_FDPS) -I$(PS_PATH_PIKG)

####################
###   Compiler   ###
####################
CC = clang++ -std=c++11
#CC = mpic++ -std=c++11
#CC = time FCCpx -Nclang
#CC = time mpiFCCpx -Nclang

########################
###   MPI & OpenMP   ###
########################
#use_OpenMP = yes
#use_MPI = yes

################
###   SIMD   ###
################
#use_AVX2 = yes
#use_AVX512 = yes
#use_SVE = yes

#########################
###   Other Options   ###
#########################
individual_cutoff = yes
quad = yes
#isotropic = yes
#gas_drag = yes
#integrate_6th = yes
#absorb_particle = yes
#test_particle = yes

collision = yes
#Kominami = yes
#Chembers = yes
#merge_binary = yes

#output_detail = yes
#calc_time = yes


CFLAGS = -O2
CFLAGS += -Wall
CFLAGS += -funroll-loops

ifeq ($(use_OpenMP),yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif
ifeq ($(use_MPI),yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif

CFLAGS_KERNEL_EPEP = --kernel-name CalcForceLongEPEP --epi-name EPGrav --epj-name EPGrav --force-name ForceGrav
CFLAGS_KERNEL_EPSP = --kernel-name CalcForceLongEPSP --epi-name EPGrav --epj-name SPGrav --force-name ForceGrav

ifeq ($(use_AVX2),yes)
#CFLAGS += -march=core-avx2
CFLAGS += -mavx2 -mfma -ffast-math
CFLAGS_KERNEL_EPEP += --conversion-type AVX2
CFLAGS_KERNEL_EPSP += --conversion-type AVX2
else
ifeq ($(use_AVX512),yes)
#CFLAGS += -march=skylake-avx512
CFLAGS += -mavx512f -mavx512dq -ffast-math
CFLAGS_KERNEL_EPEP += --conversion-type AVX-512
CFLAGS_KERNEL_EPSP += --conversion-type AVX-512
else
ifeq ($(use_SVE),yes)
CFLAGS += -Ofast
CFLAGS_KERNEL_EPEP += --conversion-type A64FX
CFLAGS_KERNEL_EPSP += --conversion-type A64FX
endif
endif
endif


ifeq ($(individual_cutoff),yes)

CFLAGS += -DUSE_INDIVIDUAL_CUTOFF
KERNEL_EPEP = gravity_kernel_EPEP_individual.pikg
ifeq ($(quad),yes)
CFLAGS += -DUSE_QUAD
KERNEL_EPSP = gravity_kernel_EPSP_quad_individual.pikg
#CFLAGS_KERNEL_EPSP += --epj-name PS::SPJQuadrupoleInAndOut
else
KERNEL_EPSP = gravity_kernel_EPSP_mono_individual.pikg
#CFLAGS_KERNEL_EPSP += --epj-name PS::SPJMonopoleInAndOut
endif

else

KERNEL_EPEP = gravity_kernel_EPEP_shared.pikg
ifeq ($(quad),yes)
CFLAGS += -DUSE_QUAD
KERNEL_EPSP = gravity_kernel_EPSP_quad_individual.pikg
#CFLAGS_KERNEL_EPSP += --epj-name PS::SPJQuadrupoleScatter
else
KERNEL_EPSP = gravity_kernel_EPSP_mono_shared.pikg
#CFLAGS_KERNEL_EPSP += --epj-name PS::SPJMonopoleScatter
endif

endif

ifeq ($(isotropic),yes)
CFLAGS += -DISOTROPIC
endif
ifeq ($(gas_drag),yes)
CFLAGS += -DGAS_DRAG
endif
ifeq ($(integrate_6th),yes)
CFLAGS += -DINTEGRATE_6TH_SUN
endif
ifeq ($(absorb_particle),yes)
CFLAGS += -DABSORB_PTCL
endif
ifeq ($(test_particle),yes)
CFLAGS += -DTEST_PTCL
endif

ifeq ($(collision),yes)
CFLAGS += -DCOLLISION
endif
ifeq ($(Kominami),yes)
CFLAGS += -DKOMINAMI
else
ifeq ($(Chambers),yes)
CFLAGS += -DCHAMBERS
endif
endif
ifeq ($(merge_binary),yes)
CFLAGS += -DMERGE_BINARY
endif

ifeq ($(output_detail),yes)
CFLAGS += -DOUTPUT_DETAIL
endif
ifeq ($(calc_time),yes)
CFLAGS += -DCALC_WTIME
endif

CFLAGS += -DMONAR
#CFLAGS += -pg

SRC = main_p3t.cpp
PROGRAM = gplum.out
HEADER = mathfunc.h kepler.h energy.h particle.h disk.h gravity.h gravity_kernel.hpp collisionA.h collisionB.h hermite.h hard.h read.h time.h func.h


$(PROGRAM): $(SRC) gravity_kernel_EPEP.hpp gravity_kernel_EPSP.hpp
	$(CC) -MMD $(INC) $(CFLAGS) -o $@ $(SRC)

gravity_kernel_EPEP.hpp: $(KERNEL_EPEP)
	../PIKG/bin/pikg -i $(KERNEL_EPEP) -o gravity_kernel_EPEP.hpp $(CFLAGS_KERNEL_EPEP)
	sed -i -e 's/std::cout/\/\/std::cout/g' gravity_kernel_EPEP.hpp
	sed -i -e 's/std::cerr/\/\/std::cerr/g' gravity_kernel_EPEP.hpp	

gravity_kernel_EPSP.hpp: $(KERNEL_EPSP)
	../PIKG/bin/pikg -i $(KERNEL_EPSP) -o gravity_kernel_EPSP.hpp $(CFLAGS_KERNEL_EPSP)
	sed -i -e 's/quad__DOT__/quad./g' gravity_kernel_EPSP.hpp
	sed -i -e 's/std::cout/\/\/std::cout/g' gravity_kernel_EPSP.hpp
	sed -i -e 's/std::cerr/\/\/std::cerr/g' gravity_kernel_EPSP.hpp	

all: clean $(PROGRAM)

clean:
	rm -f *.out *.o *.d *~
	rm -f gravity_kernel*.hpp*

