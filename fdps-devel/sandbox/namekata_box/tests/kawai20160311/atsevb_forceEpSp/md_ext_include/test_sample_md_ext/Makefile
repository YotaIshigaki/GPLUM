#--------------------------------------------------------------------
#  Makefile for test function of <md_ext_***.hpp>.
#    2016/2/25  K.Kawai Tokumasu-Lab.
#--------------------------------------------------------------------

#--- C++ compile target
SRC_CPP = cpp_test_multi_table.cpp
OBJ_CPP = $(SRC_CPP:%.cpp=%.o)
OBJS    = $(OBJ_CPP) $(OBJ_F)
PROGRAM = test.out

#--- PATH for FDPS src directry
#------ for paolo in Tokumasu Lab.
PS_DIR   = /home/kawai/local/FDPS_2.0_20160125
#PS_DIR   = /home/kawai/local/FDPS_2.0_20160217
#PS_DIR   = /home/kawai/local/FDPS_1.1_20150828

PS_PATH  = -I$(PS_DIR)/src/

#--- PATH for <md_ext.hpp>
MD_EXT_PATH = -I../md_ext_include/

#--- compiler settings
CC = g++
#CC = clang
#CC = mpicxx
#CFLAGS = -O3 -ffast-math -funroll-loops -lm -std=c++11 -g3
CFLAGS = -O0 -lm -std=c++11 -g3 -Wall

#--- enable debug flags for <md_ext.hpp>
CFLAGS += -DCOEF_TABLE_WINDOW_TEST

#--- dependency check
DEPS = $(SRC_CPP:%.cpp=%.d)


all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(PS_PATH) $(MD_EXT_PATH) $(LIB_FFTW) $(CFLAGS) -o $@ $< $(LIB_PM)

$(OBJ_CPP): $(SRC_CPP)
	$(CC) -c $(PS_PATH) $(INCLUDE_FFTW) $(CFLAGS) -MMD $<

    
.PHONY: clean
clean:
	rm $(PROGRAM) $(DEPS) $(OBJS)


-include $(DEPS)
    
