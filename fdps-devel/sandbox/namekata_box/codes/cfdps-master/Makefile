#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
#FDPS_LOC = ../../../
FDPS_LOC = ../../src/fdps/fdps/
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_FTN_MOD_DIR = $(FDPS_LOC)/src/fortran_interface/modules
FDPS_FTN_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_ftn_if.py
EXPORTSRCS = Makefile\
             user_defined.c\
             user_defined.h\
             FDPS_vector.h\
             README.md\
             convert_c_struct_to_f90.rb\
             convert_f90_struct_to_c.rb\
             convert_f90_if_to_c.rb\
             main.c
EXPORTDIR = export

# (ii) Variables to specify compilers and compile options
# Serial or OpenMP cases
FC=gfortran
CXX=g++
CC=gcc
# MPI case
#FC=mpif90
CXX=mpic++
# [Option 1] w/o optimization
#FCFLAGS = -std=f2003 -O0 -Wall
CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
FCFLAGS = -std=f2003 -O3 -ffast-math -funroll-loops -finline-functions
CXXFLAGS = -O3 -ffast-math -funroll-loops $(FDPS_INC)
LDFLAGS = -lgfortran 
# OpenMP options
#FCFLAGS  += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# MPI options
#FCFLAGS  += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

#-----------------------------------------------------------------------
#   Source files
#-----------------------------------------------------------------------
%.o : %.F90
	$(FC) $(FCFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

SRC_USER_DEFINED_TYPE = user_defined.F90
SRC_USER = f_main.F90
SRC_FDPS_MOD = $(wildcard $(FDPS_FTN_MOD_DIR)/*.F90)
SRC_FTN = $(SRC_FDPS_MOD) \
	  $(SRC_USER_DEFINED_TYPE) \
	  FDPS_module.F90 \
	  $(SRC_USER)
SRC_CXX = FDPS_ftn_if.cpp \
	  FDPS_Manipulators.cpp \
	  mainf.cpp

OBJ_USER_DEFINED_TYPE	= $(SRC_USER_DEFINED_TYPE:F90=o)
OBJ_USER		= $(SRC_USER:F90=o)
OBJ_FDPS_MOD		= $(notdir $(SRC_FDPS_MOD:F90=o))
OBJ_FTN			= $(notdir $(SRC_FTN:F90=o))
OBJ_CXX			= $(SRC_CXX:cpp=o)
OBJ			= $(OBJ_FTN) $(OBJ_CXX)

VPATH = $(FDPS_FTN_MOD_DIR)
TARGET = nbody.out

$(TARGET): $(OBJ) result
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)

result:
	mkdir result

user_defined.F90: user_defined.h Makefile
	ruby convert_c_struct_to_f90.rb  user_defined.h> user_defined.F90
$(SRC_CXX) FDPS_module.F90: $(SRC_USER_DEFINED_TYPE)
	$(FDPS_FTN_IF_GENERATOR) $(SRC_USER_DEFINED_TYPE) --output ./ 

FDPS_super_particle.o: FDPS_vector.o FDPS_matrix.o

$(OBJ_USER_DEFINED_TYPE): $(OBJ_FDPS_MOD)

FDPS_module.o: $(OBJ_USER_DEFINED_TYPE)

$(OBJ_USER): $(OBJ_USER_DEFINED_TYPE) FDPS_module.o

clean:
	rm -f *.o *.s *.mod $(TARGET) *.dat

distclean: clean
	rm -f $(SRC_CXX) FDPS_Manipulators.h  FDPS_module.F90 user_defined.hpp 
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)

user_defined.o : user_defined.c
	$(CC) $(CXXFLAGS) -c user_defined.c
FDPS_c_if.h:   FDPS_ftn_if.cpp  convert_f90_if_to_c.rb
	ruby convert_f90_if_to_c.rb FDPS_ftn_if.cpp > FDPS_c_if.h
FDPS_types.h:   $(FDPS_FTN_MOD_DIR)/FDPS_time_profile.F90 $(FDPS_FTN_MOD_DIR)/FDPS_super_particle.F90 convert_f90_struct_to_c.rb
	cpp -E $(FDPS_FTN_MOD_DIR)/FDPS_super_particle.F90 .ftmp.F90
	cat .ftmp.F90  $(FDPS_FTN_MOD_DIR)/FDPS_time_profile.F90| ruby  convert_f90_struct_to_c.rb> FDPS_types.h
main.o: main.c user_defined.c user_defined.h FDPS_vector.h FDPS_c_if.h FDPS_types.h Makefile
	$(CC) $(CXXFLAGS) -c main.c
cmain.o : cmain.cpp
COBJS =   FDPS_ftn_if.o FDPS_Manipulators.o main.o  cmain.o user_defined.o
CLIBS =
fdpsc:  $(COBJS) $(CLIBS) Makefile
	$(CXX) $(COBJS)  -O3 -ffast-math -funroll-loops -I../fdps/fdps//src  -o fdpsc -L. 

exports : $(EXPORTSRCS)
	cp -p $(EXPORTSRCS) $(EXPORTDIR)

git :
	git add --all
	git commit
	git push origin master
