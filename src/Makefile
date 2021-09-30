PS_PATH_FDPS = ../FDPS/src/
PS_PATH_PIKG = ../PIKG/
INC = -I$(PS_PATH_FDPS) -I$(PS_PATH_PIKG)inc

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
quadrupole = yes
#polar_coordinate = yes
#isotropic = yes
#gas_drag = yes
#integrate_6th = yes
#absorb_particle = yes
#test_particle = yes

collision = yes
#Kominami = yes
#Chembers = yes
#Shibata = yes
#merge_binary = yes

#without_sun = yes
#constant_random_velocity = yes
#correct_interaction_gravity = yes

output_detail = yes
#calc_time = yes

p2p_fast = yes


CFLAGS = -O2
CFLAGS += -Wall
CFLAGS += -funroll-loops

CFLAGS += -DFOR_PIKG01
CFLAGS += -DPARTICLE_SIMULATOR_USE_STD_SORT

ifeq ($(use_OpenMP),yes)
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif
ifeq ($(use_MPI),yes)
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
endif

CFLAGS_KERNEL_EPEP = --kernel-name CalcForceLongEPEP --epi-name EPI_t --epj-name EPJ_t --force-name Force_t
CFLAGS_KERNEL_EPSP = --kernel-name CalcForceLongEPSP --epi-name EPI_t --epj-name SPJ_t --force-name Force_t

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
endif
ifeq ($(quadrupole),yes)
CFLAGS += -DUSE_QUAD
endif
ifeq ($(polar_coordinate),yes)
CFLAGS += -DUSE_POLAR_COORDINATE
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
else
ifeq ($(Shibata),yes)
CFLAGS += -DSHIBATA
endif
endif
endif
ifeq ($(merge_binary),yes)
CFLAGS += -DMERGE_BINARY
endif
ifeq ($(without_sun),yes)
CFLAGS += -DWITHOUT_SUN
endif
ifeq ($(constant_random_velocity),yes)
CFLAGS += -DCONSTANT_RANDOM_VELOCITY
endif
ifeq ($(correct_interaction_gravity),yes)
CFLAGS += -DCORRECT_INTERACTION_GRAVITY
endif

ifeq ($(output_detail),yes)
CFLAGS += -DOUTPUT_DETAIL
endif
ifeq ($(calc_time),yes)
CFLAGS += -DCALC_WTIME
endif

ifeq ($(p2p_fast),yes)
CFLAGS += -DUSE_P2P_FAST
endif

CFLAGS += -DMONAR
#CFLAGS += -pg

SRC = main_p3t.cpp
PROGRAM = gplum.out
HEADER = mathfunc.h kepler.h energy.h particle.h disk.h gravity.h gravity_kernel.hpp collisionA.h collisionB.h hermite.h hard.h read.h time.h func.h


$(PROGRAM): $(SRC) gravity_kernel_epep.hpp gravity_kernel_epsp.hpp
	$(CC) -MMD $(INC) $(CFLAGS) -o $@ $(SRC)

gravity_kernel_epep.hpp: gravity_kernel_epep.pikg
	./removeif.sh gravity_kernel_epep.pikg gravity_kernel_epep_noif.pikg $(CFLAGS)
	$(PS_PATH_PIKG)bin/pikg -i gravity_kernel_epep_noif.pikg -o gravity_kernel_epep.hpp $(CFLAGS_KERNEL_EPEP)
#	sed -i -e 's/std::cout/\/\/std::cout/g' gravity_kernel_epep.hpp
#	sed -i -e 's/std::cerr/\/\/std::cerr/g' gravity_kernel_epep.hpp	

gravity_kernel_epsp.hpp: gravity_kernel_epsp.pikg
	./removeif.sh gravity_kernel_epsp.pikg gravity_kernel_epsp_noif.pikg $(CFLAGS)
	$(PS_PATH_PIKG)bin/pikg -i gravity_kernel_epsp_noif.pikg -o gravity_kernel_epsp.hpp $(CFLAGS_KERNEL_EPSP)
	sed -i -e 's/quad__DOT__/quad./g' gravity_kernel_epsp.hpp
#	sed -i -e 's/std::cout/\/\/std::cout/g' gravity_kernel_epsp.hpp
#	sed -i -e 's/std::cerr/\/\/std::cerr/g' gravity_kernel_epsp.hpp	

all: clean $(PROGRAM)

clean:
	rm -f *.out *.o *.d *~
	rm -f gravity_kernel*.hpp* gravity_kernel*_noif.pikg*

