#-----------------------------------------------------------------------
#   Compiler and complile options                                       
#-----------------------------------------------------------------------
#HOST = mac
HOST = armagnac0

# [1] Mac Env.
ifeq ($(HOST),mac)
FC=gfortran-mp-5
CC=gcc
CXX=g++-mp-5

FDPS_LOC = /Users/name/fdps
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_LIB = 

FCFLAGS = -O0 -Wall
CFLAGS =
CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
LDFLAGS = -lgfortran $(FDPS_LIB)
endif

# [2] armagnac0
ifeq ($(HOST),armagnac0)
FC=gfortran
CC=gcc
CXX=g++

FDPS_LOC = /nfshome/namekata/codes/fdps
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_LIB = 

FCFLAGS = -std=f2003 -O0 -Wall
CFLAGS =
CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
LDFLAGS = -lgfortran $(FDPS_LIB)
endif

#-----------------------------------------------------------------------
#   Source files                                                        
#-----------------------------------------------------------------------
SRC_FTN = f_main.F90
SRC_C   =
SRC_CXX = main.cpp
SRC     = $(SRC_FTN) $(SRC_C) $(SRC_CXX)

OBJ_FTN = $(addsuffix .o, $(basename $(SRC_FTN))) # $(SRC1:.F90=.o)
OBJ_C   = $(addsuffix .o, $(basename $(SRC_C)))   # $(SRC1:.c=.o)
OBJ_CXX = $(addsuffix .o, $(basename $(SRC_CXX))) # $(SRC3:.cpp=.o)
OBJ     = $(OBJ_FTN) $(OBJ_C) $(OBJ_CXX)

TARGET = a.out

$(TARGET): Makefile $(SRC)
	make obj

obj: Makefile $(SRC)
	make obj_ftn
	#make obj_c
	make obj_cxx
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)
	strip $(TARGET)


obj_ftn: $(SRC_FTN) Makefile
	$(FC) $(FCFLAGS) -c $(SRC_FTN)

obj_c: $(SRC_C) Makefile
	$(CC) $(CFLAGS) -c $(SRC_C)

obj_cxx: $(SRC_CXX) Makefile
	$(CXX) $(CXXFLAGS) -c $(SRC_CXX)

clean:
	rm -f *.o *.mod $(TARGET) *.txt

