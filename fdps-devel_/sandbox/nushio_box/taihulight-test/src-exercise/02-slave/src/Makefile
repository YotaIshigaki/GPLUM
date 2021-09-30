CC = sw5cc.new
all: main.out

CFLAGS = -O3 -msimd
LIBFLAGS = -hybrid #-lslave -allshare

main.out:master.o slave.o
	$(CC) $(LIBFLAGS) master.o slave.o -o main.out

master.o:master.c
	$(CC) $(CFLAGS) -host -E master.c > master.e
	#	$(CC) $(CFLAGS) -host -s master.c -o master.s
	$(CC) $(CFLAGS) -host -c master.c -o master.o

slave.o:slave.c
	$(CC) $(CFLAGS) -slave -E slave.c > slave.e
	#$(CC) $(CFLAGS) -slave -s slave.c -o slave.s
	$(CC) $(CFLAGS) -slave -c slave.c -o slave.o

run:
	bsub -I -b -q q_sw_expr -n 1 -cgsp 64 -share_size 4096 -host_stack 128 ./main.out

core_functions.o: core_functions.S
	$(CC) $(CFLAGS) -slave -c $^ -o $@
