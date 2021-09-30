#==========================
# Calculation conditions
#==========================
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
CC = sw5cc.new
all: main.out

CFLAGS = -O2 -Wall -msimd
LIBFLAGS = -hybrid -lm #-lslave -allshare

TARGET = prefix_sum.x
MY_DIR = /home/export/online1/riken/sandbox/namekata_box/codes/taihulight-jobs/prefix_sum

$(TARGET): master.o slave.o
	$(CC) master.o slave.o $(LIBFLAGS) -o $(TARGET)

master.o: master.c
	$(CC) $(CFLAGS) -host -E master.c > master.e
	$(CC) $(CFLAGS) -host -S master.c -o master.s
	$(CC) $(CFLAGS) -host -c master.s -o master.o

slave.o: slave.s
	$(CC) $(CFLAGS) -slave -c slave.s -o slave.o

slave.s: slave.c
#	$(CC) $(CFLAGS) -slave -E slave.c > slave.e
	$(CC) $(CFLAGS) -slave -S slave.c -o slave.s

clean:
	rm -f $(TARGET) *.o *.s *.e

ssh_setting:
	ssh-keygen -R $(ACCESS_POINT)

send_src: ssh_setting
	$(SSH) $(ACCESS_POINT) "mkdir -p $(MY_DIR)"
	$(RSYNC_SEND) --include='Makefile' --include="*.h" --include="*.c" --exclude='*' ./ $(ACCESS_POINT):$(MY_DIR)

compile: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); make clean; make -B $(TARGET)"

submit: ssh_setting
	$(SSH) $(ACCESS_POINT) "cd $(MY_DIR); \
	bsub -I -b -q q_sw_expr -n 1 -cgsp 64 -share_size 6500 -host_stack 512 ./$(TARGET) > $(LOGFILE) 2>&1"

get: ssh_setting
	$(RSYNC_RECV) $(ACCESS_POINT):$(MY_DIR)/00stdout.log ./

scs: send_src compile submit

run: main.out
	bsub -I -b -q q_sw_expr -n 1 -cgsp 64 -share_size 4096 -host_stack 128 ./main.out

