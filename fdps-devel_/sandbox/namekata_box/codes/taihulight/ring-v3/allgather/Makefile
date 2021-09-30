#==========================
# Calculation conditions
#==========================
NRAD = 1
NPHI = 8
NPROCS = $(shell expr $(NRAD) \* $(NPHI))
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
CC = sw5cc
CXX = swg++ -I/usr/sw-mpp/mpi2/include -I/usr/sw-mpp/swcc/sw5gcc-binary/include
CLINK = mpiCC -hybrid

TARGET = allgather-namekata.out
MY_DIR = /home/export/online1/riken/sandbox/namekata_box/codes/taihulight-jobs/allgather

CFLAGS = -O3
#CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
#CFLAGS += -fopenmp

OBJS = diag_gather.o 

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CLINK) $(INC) $(CFLAGS) $^ $(CLIBS) -o $@

diag_gather.o: diag_gather.cpp 
	$(CXX) $(INC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f $(TARGET) *.o *~ 

ssh_setting:
	ssh-keygen -R $(ACCESS_POINT)

send_src: ssh_setting
	$(RSYNC_SEND) --include='Makefile' --include="*.cpp" --exclude='*' ./ $(ACCESS_POINT):$(MY_DIR)

compile: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); make clean; make -B $(TARGET)"

submit: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q $(QUEUE) -n $(NPROCS) -cgsp 64 -share_size 6500 -host_stack 512 ./$(TARGET) $(NRAD) $(NPHI) > $(LOGFILE) 2>&1"

get: ssh_setting
	$(RSYNC_RECV) $(ACCESS_POINT):$(MY_DIR)/00stdout.log ./
	$(RSYNC_RECV) $(ACCESS_POINT):$(MY_DIR)/*.txt ./

sc: send_src compile

scs: send_src compile submit

