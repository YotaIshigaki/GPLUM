HOST = $(shell hostname)

SRC_FDPS_DIR = /home/masaki/project/fdps.git/src/
CXX = mpicxx
CXXFLAGS = -O3
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -Wall
CXXFLAGS += -std=c++1z
CXXFLAGS += -DDEV_CODE
CXXFLAGS += -DADD_RNADOM_SHIFT
#CXXFLAGS += -DEXPAND_BOX

DEST_SYS = hokusai

ifeq ($(findstring hokusai, $(DEST_SYS)), hokusai)
USER = iwasawa
DEST_DIR = /home/iwasawa/work/paper/ring/
DEST_FDPS_DIR = $(DEST_DIR)/fdps/src
endif

ifeq ($(findstring hokusai, $(HOST)), hokusai)
CXX = mpiicpc
#DEST_DIR = /home/iwasawa/work/paper/ring/
#DEST_FDPS_DIR = $(DEST_DIR)/fdps/src
FDPS_PATH = /home/iwasawa/work/paper/ring/fdps/src
PHANTOM_PATH = /home/iwasawa/project/git/phantom-grape/PG5/newton/libpg5_avx512/
CXXFLAGS += -xCORE-AVX512
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -qopenmp
#CXXFLAGS += -xCOMMON-AVX512
#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY=1
#CXXFLAGS += -Qvec-report2
else ifeq ($(findstring DESKTOP, $(HOST)), DESKTOP)
CXX = mpicxx
#DEST_FDPS_DIR = $(HOME)/project/fdps.git/src/
FDPS_PATH = $(HOME)/project/fdps.git/src/
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
endif

INC = -I$(FDPS_PATH)
INC += -I$(PHANTOM_PATH)
LIB = -L


main_grav_mono_phantom.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DGRAV_ONLY -DENABLE_PHANTOM_GRAPE_X86 main.cpp -o $@

main_rotate_mono.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DROTATE main.cpp -o $@

main_grav_rotate_mono.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DGRAV_ONLY -DROTATE main.cpp -o $@

main_mono.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) main.cpp -o $@

main_mono_expand.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DEXPAND_BOX main.cpp -o $@

main_quad.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DQUAD main.cpp -o $@

main_grav_mono.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DGRAV_ONLY main.cpp -o $@

main_grav_quad.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DGRAV_ONLY -DQUAD main.cpp -o $@

main_grav_mono_expand.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DEXPAND_BOX -DGRAV_ONLY main.cpp -o $@

main_grav_quad_expand.out: main.cpp
	$(CXX) $(INC) $(CXXFLAGS) -DEXPAND_BOX -DGRAV_ONLY -DQUAD main.cpp -o $@

main.s: main.cpp
	$(CXX) -S $(INC) $(CXXFLAGS) main.cpp -o main.s

send_src:
	ssh $(USER)@$(DEST_SYS) "cd ./work/paper/ring; mkdir -pv fdps; cd fdps; mkdir -pv src;"
	rsync -av --include='*.hpp' --include='*.cpp' --exclude=$(SRC_FDPS_DIR)/* $(SRC_FDPS_DIR)/ $(USER)@$(DEST_SYS):$(DEST_FDPS_DIR)
	rsync -av --include='*.hpp' --include='*.cpp' --include='Makefile' \
	--exclude='.git' --exclude='*~' --exclude='*.eps' --exclude='data/' --exclude='tmp/' --exclude='result/' ../ring02 $(USER)@$(DEST_SYS):$(DEST_DIR)
