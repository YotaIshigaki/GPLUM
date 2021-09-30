#====================================================================
#  Makefile for "md_fdps"
#====================================================================

#--- compiler settings
#CXX = mpicxx
CXX = mpicxx-openmpi-gcc6

#CPPFLAGS = -lm -std=c++11 -MMD -lefence
CPPFLAGS = -lm -std=c++11 -MMD
#CPPFLAGS += -O3 -ffast-math -funroll-loops
CPPFLAGS += -O0 -Wall -Wextra -g3

#--- parallelization flag for FDPS
#CPPFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CPPFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CPPFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY=1

#--- C++ compile target
SRC_DIR = ./src
SRC_CPP = $(SRC_DIR)/md_init.cpp $(SRC_DIR)/md_fdps.cpp
INCLUDE = -I$(SRC_DIR)/ -I./generic_ext/

#--- source by automatic code generator
SRC_AUTO = $(SRC_DIR)/enum_model.hpp


#--- PATH for FDPS library
#------ default
PS_DIR = $(FDPS_ROOT)
#------ absolute path
#PS_DIR = $(HOME)/FDPS/FDPS_4.0a_20171117
#PS_DIR = $(HOME)/FDPS/FDPS_4.0a_gcc54

#--- PATH for FFTw library
#------ default
FFTW_DIR = $(FFTW_ROOT)
#------ absolute path
#FFTW_DIR = $(HOME)/local/fftw-3.3.6
#FFTW_DIR = $(HOME)/local/fftw-3.3.6-ompi-2.1.1-gcc-5.4

#--- PATH for Google test framework
GTEST_DIR = $(GTEST_ROOT)


#--- preparing for PATH
PS_PATH  = -I$(PS_DIR)/src/
PS_PATH += -I$(PS_DIR)/src/particle_mesh/
LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm.a
#LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm_debug.a

INCLUDE_FFTW  = -I$(FFTW_DIR)/include/
LIB_FFTW      = -L$(FFTW_DIR)/lib/ -lfftw3f_mpi -lfftw3f
LIB_FFTW_STATIC  = $(FFTW_DIR)/lib/libfftw3f_mpi.a
LIB_FFTW_STATIC += $(FFTW_DIR)/lib/libfftw3f.a

INCLUDE += $(PS_PATH) $(INCLUDE_FFTW)

#LIBS = $(LIB_PM) $(LIB_FFTW)
LIBS = $(LIB_PM) $(LIB_FFTW_STATIC)


#--- preparing for main code
OBJ_DIR = obj
ifeq "$(strip $(OBJ_DIR))" ""
  OBJ_DIR = .
endif

OBJS  = $(addprefix $(OBJ_DIR)/, $(notdir $(SRC_CPP:.cpp=.o) ) )
DEPS  = $(OBJS:.o=.d)
TGT   = $(notdir $(SRC_CPP:.cpp=.x) )


#--- preparing for unit test
GTEST_SRCDIR = unit_test
GTEST_SRCS :=
REL := $(GTEST_SRCDIR)
include $(REL)/Makefile

GTEST_OBJDIR = test_obj
ifeq "$(strip $(GTEST_OBJDIR))" ""
  GTEST_OBJDIR = .
endif

GTEST_EXEDIR = test_bin
ifeq "$(strip $(GTEST_EXEDIR))" ""
  GTEST_EXEDIR = .
endif

GTEST_OBJS = $(addprefix $(GTEST_OBJDIR)/, $(notdir $(GTEST_SRCS:.cpp=.o) ) )
GTEST_DEPS = $(GTEST_OBJS:.o=.d)
GTEST_TGT  = $(addprefix $(GTEST_EXEDIR)/, $(notdir $(GTEST_SRCS:.cpp=) ) )

GTEST_INCLUDE  = $(INCLUDE)
GTEST_INCLUDE += -I$(GTEST_ROOT)/googletest/include

GTEST_LIBS  = $(LIBS)
GTEST_LIBS += $(GTEST_ROOT)/build/googlemock/gtest/libgtest.a
#GTEST_LIBS += $(GTEST_ROOT)/build/googlemock/gtest/libgtest_main.a

GTEST_FLAGS = $(CPPFLAGS) -isystem -pthread


#=======================
#  main MD program
#=======================
main: $(TGT)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_AUTO)
	@[ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(INCLUDE) -o $@ -c $<

%.x: $(OBJ_DIR)/%.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(LIBS)


#=======================
#  tests for setting / parameter file
#=======================
test_condition: ./unit_test/test_loading_condition.cpp $(SRC_AUTO)
	@[ -d $(GTEST_EXEDIR) ] || mkdir -p $(GTEST_EXEDIR)
	$(CXX) $(INCLUDE) $(CPPFLAGS) $< -o $(GTEST_EXEDIR)/test_condition
	mpirun -n 2 $(GTEST_EXEDIR)/test_condition

test_model: ./unit_test/test_loading_model.cpp $(SRC_AUTO)
	@[ -d $(GTEST_EXEDIR) ] || mkdir -p $(GTEST_EXEDIR)
	$(CXX) $(INCLUDE) $(CPPFLAGS) $< -DTEST_MOL_INSTALL -o $(GTEST_EXEDIR)/test_model
	$(GTEST_EXEDIR)/test_model  AA_wat_SPC_Fw

test_param: ./unit_test/test_param_file.cpp $(SRC_AUTO)
	@[ -d $(GTEST_EXEDIR) ] || mkdir -p $(GTEST_EXEDIR)
	$(CXX) $(INCLUDE) $(CPPFLAGS) $< -o $(GTEST_EXEDIR)/test_param
	$(GTEST_EXEDIR)/test_param  AA_C6H5_CH3


#=======================
#  unit tests
#=======================
gtest: $(GTEST_TGT)

$(GTEST_OBJDIR)/%.o: $(GTEST_SRCDIR)/%.cpp $(SRC_AUTO)
	@[ -d $(GTEST_OBJDIR) ] || mkdir -p $(GTEST_OBJDIR)
	$(CXX) $(GTEST_FLAGS) $(GTEST_INCLUDE) -o $@ -c $<

$(GTEST_EXEDIR)/%: $(GTEST_OBJDIR)/%.o
	@[ -d $(GTEST_EXEDIR) ] || mkdir -p $(GTEST_EXEDIR)
	$(CXX) $(GTEST_FLAGS) -o $@ $<  $(GTEST_LIBS)


#=======================
#  script call
#=======================
$(SRC_DIR)/enum_model.hpp:
	./script/convert_model_indicator.py

.PHONY: clean
clean:
	rm -rf $(TGT) $(OBJ_DIR)
	rm -rf $(GTEST_EXEDIR) $(GTEST_OBJDIR)

.PHONY: clean_all
clean_all:
	make clean
	rm -rf ./posdata ./pdb ./resume *.dat $(SRC_AUTO)

-include $(DEPS) $(GTEST_DEPS)
