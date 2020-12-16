NSCCWX = riken@42.0.1.125
PS_PATH = ./src/
CPE_KERNEL = ./cpe_kernel/
INC = -I$(PS_PATH) -I$(CPE_KERNEL)

#CC  = time gcc
#CXX = mpicxx
#CLINK=mpicxx

CC = sw5cc
CXX = swg++ -I/usr/sw-mpp/mpi2/include -I/usr/sw-mpp/swcc/sw5gcc-binary/include
CLINK=mpiCC -hybrid

CFLAGS = -O3
#CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CFLAGS += -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CFLAGS += -DPARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
#CFLAGS += -DPARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
#CFLAGS += -DCOLLISION_TEST
CFLAGS += -DCANCEL_FORCE
CFLAGS += -DUSE_STD_SORT
CFLAGS += -DROTATE_PTCLS
#CFLAGS += -DCHECK_NAN
CFLAGS += -DSUNWAY
CFLAGS += -DSUNWAY_FORCE_KERNEL
#CFLAGS += -DMPE_FORCE_KERNEL
CFLAGS += -DSPLIT
#CFLAGS += -DPROFILE
#CFLAGS += -DDEBUG
#CFLAGS += -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
CFLAGS += -DSUNWAY_DEVELOP


OBJS = nbody.o copy.o gen_morton_key.o calc_moment.o rotate.o kernel_sunway_dyn.o force_loop.o posvel_loop.o rsqrt_loop.o force_counter_loop.o spj_loop.o mpe_force_kernel.o

TARGET = nbody-namekata.out

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CLINK) $(INC) $(CFLAGS) $^ $(CLIBS) -o $@

nbody_host.out: nbody.cpp
	$(CXX) $(INC) $(CFLAGS) $^ $(CLIBS) -o $@

nbody.o: nbody.cpp 
	$(CXX) $(INC) $(CFLAGS) -c $^ -o $@

copy.o: ./cpe_kernel/copy.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

gen_morton_key.o: ./cpe_kernel/gen_morton_key.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

calc_moment.o: ./cpe_kernel/calc_moment.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

rotate.o: ./cpe_kernel/rotate.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

kernel_sunway_dyn.o: ./cpe_kernel/kernel_sunway_dyn.c
	#$(CC) -slave -msimd $(INC) $(CFLAGS) -S $^ -o kernel_sunway_dyn.s
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

force_loop.o: ./cpe_kernel/force_loop.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

posvel_loop.o: ./cpe_kernel/posvel_loop.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

rsqrt_loop.o: ./cpe_kernel/rsqrt_loop.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

force_counter_loop.o: ./cpe_kernel/force_counter_loop.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

spj_loop.o: ./cpe_kernel/spj_loop.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

mpe_force_kernel.o: mpe_force_kernel.cpp
	$(CXX) $(INC) $(CFLAGS) -c $^ -o $@

nbody.s: nbody.cpp 
	$(CXX) $(INC) $(CFLAGS) -S $^ -o $@

clean:
	rm -f $(TARGET) *.o *~ 

send_src:
	ssh $(NSCCWX) "cd ~/iwasawa_box/submit_job/; mkdir -p ./tree/; cd ./tree/; mkdir -p ./src/"
	rsync -av --include='src' --include='src/*.hpp' --include='cpe_kernel' --include='cpe_kernel/*.c' --include='cpe_kernel/*.h' --include='Makefile' --include="*.h" --include="*.hpp" --include="*.c" --include="*.cpp" --exclude='*' ./ $(NSCCWX):~/online1/sandbox/iwasawa_box/submit_job/tree/

compile:
	ssh $(NSCCWX) "cd ~/iwasawa_box/submit_job/tree/; make distclean; make -B $(TARGET)"

#submit:
#	ssh $(NSCCWX) "cd ~/iwasawa_box/submit_job/tree/; \
#	bsub -I -b -q q_sw_gb_2 -N 1024 -np 4 -cgsp 64 -share_size 4096 ./$(TARGET) -N 2000000000 -n 512 -t 0.5"

#submit:
#	ssh $(NSCCWX) "cd ~/iwasawa_box/submit_job/tree/; \
#	bsub -I -b -q q_sw_gb_1 -N 4096 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1638400000 -n 512 -t 0.5"

submit:
	ssh $(NSCCWX) "cd ~/iwasawa_box/submit_job/tree/; \
	bsub -b -q q_sw_expr -N 1 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1638400 -n 512 -t 0.5"

#submit:
#	ssh $(NSCCWX) "cd ~/iwasawa_box/submit_job/tree/; \
#	bsub -I -b -q q_sw_gb_1 -N 1024 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 20480000 -n 512 -t 0.5"

#submit:
#	ssh $(NSCCWX) "cd ~/iwasawa_box/submit_job/tree/; \
#	bsub -I -b -q q_sw_gb_1 -N 4 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 160000 -n 512 -t 0.5"

#submit:
#	ssh $(NSCCWX) "cd ~/iwasawa_box/submit_job/tree/; \
#	bsub -I -b -q q_sw_gb_2 -N 4096 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1600000000 -n 512 -t 0.5"

get:
	rsync -avx $(NSCCWX):~/iwasawa_box/submit_job/tree/result/debug_00000.dat ./

sc: send_src compile

scs: send_src compile submit

run:
	#bsub -I -b -q q_sw_cn_3 -N 1 -np 1 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 131072 -n 512 -t 0.5 > 00stdout.log
	#bsub -I -b -q q_sw_cn_3 -N 1 -np 1 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 131072 -n 512 -t 0.5 
	#bsub -I -b -q q_sw_zhengji -N 32 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1048576 -n 512 -t 0.5 > 00stdout.log
	#bsub -I -b -q q_sw_expr -N 1 -np 1 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1048576 -n 512 -t 0.5 > 00stdout.log
	#bsub -I -b -q q_sw_gb_2 -N 4 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1048576 -n 512 -t 0.5
	#bsub -I -b -q q_sw_cn_2 -N 16 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1048576 -n 512 -t 0.5
	#bsub -I -b -q q_sw_cn_2 -N 4 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 3000000 -n 512 -t 0.5
	#bsub -I -b -q q_sw_cn_3 -N 128 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 131072 -n 512 -t 0.5 >00stdout.log
	#bsub -I -b -q q_sw_gb_2 -N 512 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1048576 -n 512 -t 0.5
	#bsub -I -b -q q_sw_gb_1 -N 18144 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 100000 -n 512 -t 0.5
	#bsub -I -b -q q_sw_gb_1 -n 512 -cgsp 64 -share_size 4096 ./$(TARGET) -N 2097152 -n 512 -t 0.5
	#
	bsub -I -b -q q_sw_zhengji -N 256 -np 4 -cgsp 64 -share_size 7000 -host_stack 512 ./$(TARGET) -N 1048576 -n 512 -t 0.5 > 00stdout.log
