#==========================
# Calculation conditions
#==========================
#NPROCS = 16
#NPROCS = 32
#NPROCS = 64
NPROCS = 256
#NPROCS = 1024
#NPROCS = 2048
#NPROCS = 4096
#NPROCS = 8192
#NPROCS = 16384
#NPTCL_PER_PROC = 1024
#NPTCL_PER_PROC = 2048
#NPTCL_PER_PROC = 4096
#NPTCL_PER_PROC = 131072
#NPTCL_PER_PROC = 1048576
#NPTCL_PER_PROC = 2097152
#NPTCL_PER_PROC = 3145728
#NPTCL_PER_PROC = 4194304
#NPTCL_PER_PROC = 5242880
#NPTCL_PER_PROC = 6291456
#NPTCL_PER_PROC = 7340032
#NPTCL_PER_PROC = 8388608
#NPTCL_PER_PROC = 9437184
NPTCL_PER_PROC = 10485760 
#NPTCL_PER_PROC = 11534336
#NPTCL_PER_PROC = 12582912
#NPTCL_PER_PROC = 13631488
#NPTCL_PER_PROC = 14680064
#NPTCL_PER_PROC = 15728640
#NPTCL_PER_PROC = 16777216
N_GROUP_LIMIT = 1024
THETA_GRV = 0.5
N_LEAF_LIMIT = 16
#QUEUE = q_sw_expr
#QUEUE = q_sw_share
QUEUE = q_sw_gb
LOGFILE = 00stdout.log

#==========================
# Connection information
#==========================
USER = riken
#REMOTE_HOST_NAME = 42.0.1.125
REMOTE_HOST_NAME = 41.0.0.188
ACCESS_POINT = $(USER)@$(REMOTE_HOST_NAME)
ifeq ($(REMOTE_HOST_NAME),41.0.0.188)
SSH = ssh -o 'StrictHostKeyChecking=no'
RSYNC_SEND = rsync -ave "ssh -o StrictHostKeyChecking=no"
RSYNC_RECV = rsync -avxe "ssh -o StrictHostKeyChecking=no"
else
SSH = ssh 
RSYNC_SEND = rsync -ave "ssh"
RSYNC_RECV = rsync -avxe "ssh"
endif

#==========================
# Compile configuration
#==========================
PS_PATH = ./src/
CPE_KERNEL = ./cpe_kernel/
INC = -I$(PS_PATH) -I$(CPE_KERNEL)

CC = sw5cc
CXX = swg++ -I/usr/sw-mpp/mpi2/include -I/usr/sw-mpp/swcc/sw5gcc-binary/include
CLINK=mpiCC -hybrid

WORKER_NAME = namekata
#WORKER_NAME = iwasawa

TARGET = nbody-$(WORKER_NAME).out
TARGET2 = reproduce_check-$(WORKER_NAME).out

#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_debug_sahllow
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_debug_deep
MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_gb_160k
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_gb_120k
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_gb_80k
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_gb_60k
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_gb_40k
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_gb_20k
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_gb_10k
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/submit_job_reproduce
#MY_DIR = /home/export/base/nsccwuxi_riken/riken/$(WORKER_NAME)_box/20180410/reproduce_check


CFLAGS = -O3
#CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CFLAGS += -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CFLAGS += -DPARTICLE_SIMULATOR_PSYS_DEBUG_PRINT
#CFLAGS += -DPARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
#CFLAGS += -DGEN_MORTON_KEY_SEQ_EXEC_FOR_DEBUG
#CFLAGS += -DCOLLISION_TEST
#CFLAGS += -DCANCEL_FORCE
CFLAGS += -DROTATE_PTCLS
CFLAGS += -DSUNWAY
CFLAGS += -DSUNWAY_FORCE_KERNEL
#CFLAGS += -DMPE_FORCE_KERNEL
CFLAGS += -DUSE_DMA
#CFLAGS += -DSPLIT
CFLAGS += -DFULL_ASM
#CFLAGS += -DPROFILE
#CFLAGS += -DDEBUG
CFLAGS += -DSUNWAY_DEVELOP
#CFLAGS += -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
CFLAGS += -DUSE_NITADORI_ALLGATHER_IN_EX_LET

######### DEBUG OPTION (from here) ###########

CFLAGS += -DCHECK_ANGULAR_MOMENTUM

# 以下、4つが消えている状態が測定用
CFLAGS += -DPARTICLE_SIMULATOR_IMPL3_DEBUG_PRINT
CFLAGS += -DDEBUG_PRINT_MAIN
CFLAGS += -DREPRODUCE_CHECK
CFLAGS += -DPOLLING_FOR_CPE_ERROR
# POLLING_FOR_CPE_ERRORだけが消えている状態がsemi-debugモード

####CFLAGS += -DPRINT_HOST_NAME


####CFLAGS += -DDUMP_TREE_CELL_GLOBAL
####CFLAGS += -DDUMP_MEMORY_SIZE
####CFLAGS += -DCHECK_ANGULAR_MOMENTUM_2
####CFLAGS += -DCHECK_NAN
####CFLAGS += -DCHECK_RMAX_RMIN


####CFLAGS += -DDEBUG_PRINT_DISPATCH_KERNEL
####CFLAGS += -DDEBUG_SAMPLE_SORT
####CFLAGS += -DDUMP_TARGET_PTCL
####CFLAGS += -DUSE_BARRIER_FOR_MEASURE
####CFLAGS += -DDEBUG_SHIFT_CENTER_Z
####CFLAGS += -DPLANET_FORCE_ONLY


#CFLAGS += -DDEBUG_PRINT_CALC_FORCE_UNSAFE # check cm mass and cm pos of all interaction list for each ip group
#CFLAGS += -DDEBUG_PRING_CALC_FORCE_UNSAFE # check cm mass and cm pos of all interaction list for each ip group

#CFLAGS += -DDEBUG_PRINT_EX_PTCL6
#CFLAGS += -DDEBUG_EX_PTCL7 # shift particle pos to check ex_ptcl
#CFLAGS += -DDEBUG_PRINT_EX_PTCL7
#CFLAGS += -DDEBUG_EX_PTCL8 # shift particle pos to check ex_ptcl
#CFLAGS += -DDEBUG_PRINT_EX_PTCL8
####CFLAGS += -DDEBUG_PRINT_EX_PTCL8_V2
#CFLAGS += -DDEBUG_PRINT_EX_PTCL_SW
#CFLAGS += -DDEBUG_MAKELIST
#CFLAGS += -DDEBUG_PRING_CALC_FORCE_LOOP_MERGE
#CFLAGS += -DDEBUG_PRINT_SETLET_2
#CFLAGS += -DDEBUG_EXLET_SUNWAY
#CFLAGS += -DDEBUG_PRINT_EX_LET_3
#CFLAGS += -DDEBUG_PRINT_SORT_GLB_3

#CFLAGS += -DDEBUG_PRINT_EX_LET_SUPER_DOMAIN
#CFLAGS += -DDEBUG_PRINT_EX_LET_SUPER_DOMAIN_2
#CFLAGS += -DROTATE_RANGE_CHECK
#CFLAGS += -DDEBUG_PRINT_MAKE_LIST_ID3
#CFLAGS += -DCHECK_EX_PTCL
#CFLAGS += -DINITIAL_LOAD_BALANCE
#CFLAGS += -DDEBUG_CALC_FORCE
#CFLAGS += -DDEBUG_PRINT_CALC_MOMENT_LOCAL
#CFLAGS += -DDEBUG_PRINT_EX_LET_SUPER_DOMAIN
#CFLAGS += -DDEBUG_PRINT_CALC_MOMENT_LOCAL_LEAN
#CFLAGS += -DDEBUG_PRINT_CALC_MOMENT_GLOBAL_LEAN
#CFLAGS += -DDEBUG_PRINT_CALC_MOMENT_LOCAL
#CFLAGS += -DDEBUG_PRINT_CALC_MOMENT_GLOBAL
#CFLAGS += -DDEBUG_FORCE_OUTPUT
#CFLAGS += -DBALANCE_NWALK_OUT
#CFLAGS += -DDEBUG_SORT_LOCAL
#CFLAGS += -DDEBUG_0709
#CFLAGS += -DDEBUG_PRINT_MAKE_IPG
#CFLAGS += -DDEBUG_DISPATCH_KERNEL

#CFLAGS += -DFORCE_SUM_CHECK
#CFLAGS += -DFORCE_SUM_CHECK_MPE
#CFLAGS += -DFORCE_SUM_CHECK_DIRECT

#CFLAGS += -DDUMP_SATELLITE

######### DEBUG OPTION (until here)###########

#CFLAGS += -DCANCEL_PLANET_FORCE
#CFLAGS += -DSEND_ALL_LET
CFLAGS += -DPHI_R_TREE
CFLAGS += -DUSE_SUPER_DOMAIN
CFLAGS += -DCALC_MOMENT_LEAN

####CFLAGS += -DUSE_STD_SORT
CFLAGS += -DUSE_CPE_SORT

#CFLAGS += -DNO_DOUBLE_BUFFER

CFLAGS += -DSHIFT_CENTER_Z
CFLAGS += -DKEY_3D
#CFLAGS += -DBALANCE_NWALK
CFLAGS += -DBALANCE_NWALK_ON_MPE
CFLAGS += -DUSE_MY_ATAN2
#CFLAGS += -DFAST_EX_PTCL
#CFLAGS += -DMULTI_ALLREDUCE

#CFLAGS += -DWRITE_SNP
#CFLAGS += -DENERGY_CHECK
#CFLAGS += -DFORCE_CHECK
#CFLAGS += -DINTEG_CHECK
#CFLAGS += -DENERGY_CHECK_2

CFLAGS += -DBUG_FIX

CFLAGS += -DREDUCE_MEMORY

CFLAGS += -DSTATIC_DD

CFLAGS += -DREMOVE_VERTEX

CFLAGS += -DREMOVE_EPI_SORTED

CFLAGS += -DREMOVE_EPJ_SORTED_LOC

CFLAGS += -DREMOVE_TP_LOC

CFLAGS += -DUSE_MEMORY_POOL

CFLAGS += -DUSE_96BIT_KEY

#CFLAGS += -DCHECK_AT_AFTER_FIRST_LOOP

OBJS = copy.o gen_morton_key.o link_cell.o calc_moment.o make_interact_list.o rotate.o kernel_sunway_dyn.o force_loop.o posvel_loop.o rsqrt_loop.o force_counter_loop.o spj_loop.o mpe_force_kernel.o balance_nwalk.o sample_sort.o sort.o calc_moment_lean.o epj_loop.o ex_ptcl.o my_atan.o my_atan2.o copy_msort_glb_tree_reuse.o make_ipg.o rotate_and_shift_z.o set_tp_adr_ptcl.o set_ptcl_from_buffer.o

all: $(TARGET)

$(TARGET): nbody.o $(OBJS)
	$(CLINK) $(INC) $(CFLAGS) $^ $(CLIBS) -o $@

$(TARGET2): reproduce_check.o $(OBJS)
	$(CLINK) $(INC) $(CFLAGS) $^ $(CLIBS) -o $@

nbody_host.out: nbody.cpp
	$(CXX) $(INC) $(CFLAGS) $^ $(CLIBS) -o $@

nbody.o: nbody.cpp 
	$(CXX) $(INC) $(CFLAGS) -c $^ -o $@

reproduce_check.o: reproduce_check.cpp
	$(CXX) $(INC) $(CFLAGS) -c $^ -o $@

my_atan.o: ./cpe_kernel/s_atan.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

my_atan2.o: ./cpe_kernel/e_atan2.c my_atan.o
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

copy.o: ./cpe_kernel/copy.c my_atan2.o
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

copy.s: ./cpe_kernel/copy.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -S $^ -o $@

make_ipg.o: ./cpe_kernel/make_ipg.c my_atan2.o
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

copy_msort_glb_tree_reuse.o: ./cpe_kernel/copy_msort_glb_tree_reuse.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

gen_morton_key.o: ./cpe_kernel/gen_morton_key.c my_atan2.o
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

link_cell.o: ./cpe_kernel/link_cell.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

calc_moment.o: ./cpe_kernel/calc_moment.c my_atan2.o
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

calc_moment_lean.o: ./cpe_kernel/calc_moment_lean.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

make_interact_list.o: ./cpe_kernel/make_interact_list.org.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

rotate.o: ./cpe_kernel/rotate.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

kernel_sunway_dyn.o: ./cpe_kernel/kernel_sunway_dyn.c
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

epj_loop.o: ./cpe_kernel/epj_loop.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

ex_ptcl.o: ./cpe_kernel/ex_ptcl.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

balance_nwalk.o: ./cpe_kernel/balance_nwalk.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

sample_sort.o: ./cpe_kernel/sample_sort.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

rotate_and_shift_z.o: ./cpe_kernel/rotate_and_shift_z.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

set_tp_adr_ptcl.o: ./cpe_kernel/set_tp_adr_ptcl.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

set_ptcl_from_buffer.o: ./cpe_kernel/set_ptcl_from_buffer.c
	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

sort.o: ./cpe_kernel/sort.cpp
	$(CXX) $(INC) $(CFLAGS) -c $^ -o $@

#my_atan.o: ./myAtan.c
#	$(CC) -slave -msimd $(INC) $(CFLAGS) -c $^ -o $@

mpe_force_kernel.o: mpe_force_kernel.cpp
	$(CXX) $(INC) $(CFLAGS) -c $^ -o $@



nbody.s: nbody.cpp 
	$(CXX) $(INC) $(CFLAGS) -S $^ -o $@

clean:
	rm -f $(TARGET) *.o *~ 

ssh_setting:
	ssh-keygen -R $(ACCESS_POINT)

remove_src: ssh_setting
	$(SSH) $(ACCESS_POINT) "mkdir -p $(MY_DIR); cd $(MY_DIR); rm -f cpe_kernel/*.*; rm -f result/*.*; rm -f src/*.*; rm -f ./*"

send_src: ssh_setting
	$(SSH) $(ACCESS_POINT) "mkdir -p $(MY_DIR); cd $(MY_DIR); \
	mkdir -p ./src/; mkdir -p ./cpe_kernel/; mkdir -p ./result; \
	rm -f result/*.*"
	$(RSYNC_SEND) --include='src' --include='src/*.hpp' --include='cpe_kernel' --include='cpe_kernel/*.c' --include='cpe_kernel/*.h' --include='Makefile' --include="*.h" --include="*.hpp" --include="*.c" --include="*.cpp" --include="*.sh" --exclude='*' ./ $(ACCESS_POINT):$(MY_DIR)

compile: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); make distclean; make -B $(TARGET)"

compile_check: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); make distclean; make -B $(TARGET2)"

submit_gb_160k: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n 160000 -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1"

submit_gb_120k: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n 120000 -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1"

submit_gb_80k: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n 80000 -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1"

submit_gb_60k: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n 60000 -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1"

submit_gb_40k: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n 40000 -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1"

submit_gb_20k: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n 20000 -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1"

submit_gb_10k: 
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -b -q $(QUEUE) -n 10000 -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1"

submit: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n $(NPROCS) -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET) -N $(NPTCL_PER_PROC) -n $(N_GROUP_LIMIT) -t $(THETA_GRV) -l $(N_LEAF_LIMIT) > $(LOGFILE) 2>&1"

submit_check:
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n 120000 -cgsp 64 -share_size 7000 -host_stack 450 ./$(TARGET2) -N 10000000 -n 1024 -t 0.5 -l 16 > 00stdout.log 2>&1"

get_log:
	mkdir -p result
	$(RSUNC_RECV) $(ACCESS_POINT):$(MY_DIR)/00stdout.log ./result

get_snp:
	mkdir -p result
	$(RSYNC_RECV) $(ACCESS_POINT):$(MY_DIR)/result/snp_000??_00000.dat ./result
#	$(RSYNC_RECV) $(ACCESS_POINT):$(MY_DIR)/result/snp_000??_000??.dat ./result
#	$(RSYNC_RECV) $(ACCESS_POINT):$(MY_DIR)/00stdout.log ./result

get:
	$(RSYNC_RECV) $(ACCESS_POINT):$(MY_DIR)/result/debug_00000.dat ./

get_all:
	$(RSYNC_RECV) $(ACCESS_POINT):$(MY_DIR)/result ./

sc: send_src compile

scs: send_src compile submit

scsgb: send_src compile submit_gb

scsgb10k: send_src compile submit_gb_10k

scsgb20k: send_src compile submit_gb_20k

scsgb40k: send_src compile submit_gb_40k

scsgb60k: send_src compile submit_gb_60k

scsgb80k: send_src compile submit_gb_80k

scsgb120k: send_src compile submit_gb_120k

scsgb160k: send_src compile submit_gb_160k

scsd: send_src compile submit_debug

scc: send_src compile_check

sccsc: send_src compile_check submit_check

run:
	bsub -I -b -q q_sw_bcc -n 96000 -cgsp 64 -share_size 6500 -host_stack 512 ./$(TARGET) -N 1048576 -n 2048 -t 0.5 > 00stdout.log 2>&1 &

